#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <string>
#include <limits>
#include <iomanip>
#include "tryvn.h"
#include "advance.h"
#include "helper.h"
#include "icim.h"
#include <vector>
#include <algorithm>
#include <cassert>

using namespace std;

double globerr = 0.0;
double globerrvec[3], globerrvec2[3], globerrvec3[3], globerrvec4[3]; 
int globerrmx[3][3]; // glboerrmx[r] index of max err u_rr

char globorder = 2, globintorder = 3;
// if globdist = 1, use getcim345D2udist, else getcim345D2u
// note globdist == 1 means cim 3 points may use cim 5 differencing
char globdist = 0, globdistvar = 1;
char globdirectD2 = 1; //used in getcim345D2u, use cross derivative at nbr point in the same side
char globbiasstat = 1;// when approximating mixed derivative, prefer larger status?
char globcheck = 1;// check Du DDu
char globnorm = 2; // norm for residual
char globheap = 1; // used in AMG
char globoldwrongcim4 = 0, globoldwrongcim5 = 0,
     globoldcim5 = 0, globdebug = 0;


char globheat = 1;//getheatvn, usual way to calculate vn
char globexactmotion = 0, globrk = 0, globcim1order = 1, globdtdx2 = 1;

char globexactvn = 0; //use exact vn
char globregfalweno = 1;

int globGSsmooth = 0;
int globheatsmooth = 0;
double globheatsmoothtime = 0.0;

int globtestnum = 0;

char globallcim345 = 1; //use cim345 for all 2nd order point
int globcim = 345;
// 1: cim 1 
// 2: cim 2 and cim 1
// 4: cim 3 and cim 4
// 5: cim 3 and cim 5
// 345: cim 3 and combination of cim 4 and cim 5
// 0: iim

double tollinsolve = 1e-8; // tolerance of linear solver

int globlinsolve = 1;
// 0: BICGSTAB
// 1: BICGSTAB with ILU preconditioning
// 2: AMG
// 3: GMRES
// 4: hypre amg
// 5: hypre bicgstab
int globperturb = 0;
// 0: CIM perturbation
// n > 0: pertrubation using status for maximum of n steps 
// -1: no perturbation
int globsmall = 1;
// 0: not storing info for Du
// 1: storing info for Du in local index
// 2: storing info for Du in global index
int GRIDNUM = 30;//global grid number
int eindex[3] = {-1,-1,-1};
double EPSILONP = 80.0;
double EPSILONM = 2.0;
int globtesta = 0; // a term
extern int SURFOPT;


// #define FIXBANANA // surgical fix for banana shape at
bool globwritemx = false;//write coupling matrix for analysis
bool globwriteerr = false;//write error in u and Du at each grid point
ofstream outfile_mx;
ofstream outfile_uerr;
ofstream outfile_Duerr;

// ofstream delta2d("delta2d.dat", ios::out); //AGM result or dt in Advance
// ofstream epsilon2d("epsilon2d.dat", ios::out);// error of radius and normal vector component at time t
//ofstream mu2d("mu2d.dat", ios::out);// not used
//ofstream nu2d("nu2d.dat", ios::out);//not used
//ofstream eta2d("eta2d.dat", ios::out);// in checkexact, exact u value





void output(ofstream& outfile, double ***phi, int *nx)
{
   int i, thedim = 3, tindex[thedim];

   for (i = 0; i < thedim; i++)
      tindex[i] = 0;
   while (tindex[0] <= nx[0])
   {
      outfile << evalarray(phi,tindex) << " ";
      (tindex[thedim-1])++;
      for (i = thedim-1; i > 0 && tindex[i] > nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
         if (i == thedim-1)
            outfile << endl;
      }
   }
}

void output(ofstream& outfile, double **phi, int *nx)
{
   int i, thedim = 2, tindex[thedim];

   for (i = 0; i < thedim; i++)
      tindex[i] = 0;
   while (tindex[0] <= nx[0])
   {
      outfile << evalarray(phi,tindex) << " ";
      (tindex[thedim-1])++;
      for (i = thedim-1; i > 0 && tindex[i] > nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
         if (i == thedim-1)
            outfile << endl;
      }
   }
}

void output(ofstream& outfile, double **phi, char **tube, int *nx)
{
   int i, thedim = 2, tindex[thedim];

   for (i = 0; i < thedim; i++)
      tindex[i] = 0;
   while (tindex[0] <= nx[0])
   {
      if (evalarray(tube,tindex) >= 2)
         outfile << evalarray(phi,tindex) << " ";
      else
         outfile << "NaN ";
      (tindex[thedim-1])++;
      for (i = thedim-1; i > 0 && tindex[i] > nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
         if (i == thedim-1)
            outfile << endl;
      }
   }
}

void output(ofstream& outfile, double **u, double **S, int *nx)
{
   int i, thedim = 2, tindex[thedim];

   for (i = 0; i < thedim; i++)
      tindex[i] = 0;
   while (tindex[0] <= nx[0])
   {
      if (evalarray(S,tindex) < 0.0)
         outfile << evalarray(u,tindex) << " ";
      else
         outfile << "NaN ";
      (tindex[thedim-1])++;
      for (i = thedim-1; i > 0 && tindex[i] > nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
         if (i == thedim-1)
            outfile << endl;
      }
   }
}

void outputneighbors(double ***S, int *index, GridData &grid)
{
   int i, mid = 1, sindex[grid.dim], tindex[grid.dim];

   for (i = 0; i < grid.dim; i++)
      sindex[i] = 0;
   while (sindex[0] <= 2)
   {
      for (i = 0; i < grid.dim; i++)
         tindex[i] = index[i]-mid+sindex[i];
      cout << evalarray(S,tindex) << " ";
      
      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > 2; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
         cout << endl;
      }
   }
}

#if 0
void getinit(double*** &S, PBData &pb, MarchStruct &march, TempStruct &tmp, 
             GridData &grid)
{
   int i, j, r, s;
   int deg = 3;
   int tindex[grid.dim], sindex[grid.dim+1], rindex[grid.dim], qindex[grid.dim];
   double value;
   double radius, radii[grid.dim], theta;
   double x[grid.dim], y[grid.dim];

// change grid size here, bondary [-1 1]
   grid.nx[0] = GRIDNUM; //number of cell
   grid.nx[1] = grid.nx[0];
   grid.nx[2] = grid.nx[0];
   grid.a[0] = -1.0;
   grid.a[1] = -1.0;
   grid.a[2] = -1.0;
   grid.dx[0] = 2.0*fabs(grid.a[0])/grid.nx[0];
   grid.dx[1] = 2.0*fabs(grid.a[1])/grid.nx[1];
   grid.dx[2] = 2.0*fabs(grid.a[2])/grid.nx[2];
   grid.mindx = grid.dx[0];
   for (i = 1; i < grid.dim; i++)
      if (grid.dx[i] < grid.mindx)
         grid.mindx = grid.dx[i];
   grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim);
   grid.tol = 1.0e-14;
   grid.t = 0.0;

   S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

   tmp.Nfd = 10;
   tmp.fourd = new double ***[tmp.Nfd];
   tmp.fdstatus = new char[tmp.Nfd];
   for (i = 0; i < tmp.Nfd; i++)
   {
      tmp.fourd[i] = NULL;
      tmp.fdstatus[i] = 0;
   }
   tmp.A = sparseelt2ptrmatrix(grid.nx[0],grid.nx[1],grid.nx[2]); // A is 3d array,each is a pointer to sparseelt
   if (globlinsolve == 1)
      tmp.M = sparseelt2ptrmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   tmp.b = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

// change epsilon minus and plus here
   pb.epsilonm = EPSILONM;
   pb.epsilonp = EPSILONP;
   pb.gamma0 = 0.0;
   pb.psi = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   pb.dim = grid.dim;
   pb.N = 2;
   pb.x = matrix(pb.N-1,grid.dim-1);
   pb.x[0][0] = 0.1;
   pb.x[0][1] = 0.0;
   pb.x[0][2] = 0.0;
   pb.x[1][0] = -0.1;
   pb.x[1][1] = 0.0;
   pb.x[1][2] = 0.0;
   pb.Q = new double[pb.N];
   pb.c = new double[pb.N];
   pb.epsilone = new double[pb.N];
   pb.epsilonlj = new double[pb.N];
   pb.sigmalj = new double[pb.N];
   for (i = 0; i < pb.N; i++)
   {
      pb.Q[i] = 1.0;
      pb.c[i] = 1.0;
      pb.epsilone[i] = 1.0;
      pb.epsilonlj[i] = 0.0159;
      pb.sigmalj[i] = 3.653;
   }
   pb.beta = 1.0;
   pb.rho0 = 0.0333;
   pb.LJ = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

   radius = fabs(grid.a[0])/2.0;
   radii[0] = 0.5;
   radii[1] = 0.25;
   radii[2] = 0.15;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim; r++)
         x[r] = grid.a[r]+tindex[r]*grid.dx[r];
/*
      value = 0.0;
      for (r = 1; r < grid.dim; r++)
         value += x[r]*x[r];
      value = sqrt(value)-radii[0];
      theta = atan2(x[2],x[1]);
// main one used in tests
      value = sqrt(value*value+x[0]*x[0])-(0.5*(radii[1]-radii[2])*sin(3.0*theta)+
                                           0.5*(radii[1]+radii[2]));
//      value = sqrt(value*value+x[0]*x[0])-(0.5*(radii[1]-radii[2])*sin(theta)+
//                                           0.5*(radii[1]+radii[2]));
//      value = sqrt(value*value+x[0]*x[0])-radii[1];
      setvalarray(S,tindex,value);
*/
// change initial level set function S here
      value = 0.0;
      for (r = 0; r < grid.dim; r++)
//         value += x[r]*x[r]/(radii[r]*radii[r]);
         value += x[r]*x[r];
      value = sqrt(value);
//      setvalarray(S,tindex,value-1.0);
      setvalarray(S,tindex,value-0.5); //level set function of radius 0.5

      if (globperturb == 0)
         perturbelt(S,tindex,grid.tol,pb);

      setvalarray(pb.psi,tindex,0.0);
      setvalarray(pb.LJ,tindex,LJwater(tindex,pb,grid));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   if (globperturb > 0)
      perturbstatus(S,grid.tol,10,pb,grid);

   march.dist = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   march.status = cmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   march.extend = new double ***[1];
   march.extend[0] = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   march.nval = 1;
   march.dorig = S;
   march.eorig = new double***[1];
   march.eorig[0] = pb.psi;
   march.tol = 1.0e-14;
   (march.heap).head = NULL;
   (march.heap).tail = NULL;
   (march.heap).loc = heapmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   (march.heap).dim = grid.dim;
   (march.pb).epsilonm = pb.epsilonm;
   (march.pb).epsilonp = pb.epsilonp;
   (march.pb).gamma0 = pb.gamma0;
   (march.pb).psi = pb.psi;
   (march.pb).dim = grid.dim;
   (march.pb).N = pb.N;
   (march.pb).x = pb.x;
   (march.pb).Q = pb.Q;
   (march.pb).c = pb.c;
   (march.pb).epsilone = pb.epsilone;
   (march.pb).beta = pb.beta;
}
#endif




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
      project(tangent,normal,tangent,grid.dim); //tangent is projection of unit vector crossing interface to normal
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
            // if(abs(normal[r])<tol) {// if tangent = 0, that means normal is parallel to tindex2-index1, make no differece as long as tk ek = 0 in cim2
            //   tangent[r] = 1.0;
            //   break;
            // }
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

void interiorpt2(SparseElt2**** &A, double ***b, int *index, double ***S, PBData &pb, 
                 GridData &grid)
{
   int r;

   for (r = 0; r < grid.dim; r++)
      interiordirection2(A,b,index,r,S,pb,grid);
}

void interiordirection2(SparseElt2**** &A, double ***b, int *index, int r, double ***S, 
                        PBData &pb, GridData &grid)
{
   int s, j, row, col, rindex[grid.dim];
   double epsilon, value, x[grid.dim];

   if (evalarray(S,index) < 0.0)
      epsilon = pb.epsilonm;
   else
      epsilon = pb.epsilonp;
   for (s = 0; s < grid.dim; s++)
      rindex[s] = index[s];
   for (s = -1; s <= 1; s += 2)
   {
      rindex[r] = index[r]+s;
      value = -epsilon/(grid.dx[r]*grid.dx[r]);
      if (rindex[r] >= 0 && rindex[r] <= grid.nx[r])
         sparse2(index,rindex,A,value,grid);
      else
      {
         for (j = 0; j < grid.dim; j++)
            x[j] = grid.a[j]+rindex[j]*grid.dx[j];
         setvalarray(b,index,evalarray(b,index)-value*DBC(x,grid.dim,0.0));
      }
      sparse2(index,index,A,-value,grid);
   }
}
// if interior point, do nothing, leave A NULL
void interiorptsmall(SparseElt2**** &A, double ***b, int *index, double ***S, 
                     PBData &pb, GridData &grid)
{
   int r, s, rindex[grid.dim];
   double ehere, value, x[grid.dim];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   if (evalarray(S,index) < 0.0)
      ehere = pb.epsilonm;
   else
      ehere = pb.epsilonp;

   for (r = 0; r < grid.dim; r++)
   {
      for (s = -1; s <= 1; s += 2)
      {
         rindex[r] = index[r]+s;
         if (rindex[r] < 0 || rindex[r] > grid.nx[r])
          // if index is boundary. This never happens. interiorptsmall is called only at interior point
         {
            value = -ehere/(grid.dx[r]*grid.dx[r]);
            sub2coord(x,rindex,grid);
            setvalarray(b,index,evalarray(b,index)-value*DBC(x,grid.dim,0.0));
         }
      }
      rindex[r] = index[r];
   }
}

void cim1again2(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
                PBData &pb, GridData &grid)
{
   int r, s, j, k;
   int rindex[grid.dim];
   double c[grid.dim], **f; 
   double **LU, **M, value; 
   int PLR[grid.dim], PLC[grid.dim];
   double alpha, beta, tangent[grid.dim], normal[grid.dim], temp[grid.dim];
   double ethere, ehere, ebar;
// ADDING
   double sigma, d[grid.dim], thesign, tau, Dtau[grid.dim];

   LU = matrix(grid.dim-1,grid.dim-1);
   M = matrix(grid.dim-1,grid.dim-1);
   f = matrix(grid.dim-1,grid.dim-1);

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
            f[r][j] = 0.0;
         c[r] = 0.0;
         d[r] = 0.0;
      }
      for (r = 0; r < grid.dim; r++)
      {
         rindex[r] = index[r]+s;
//         if (index[0] == 24 && index[1] == 28 && index[2] == 39)
//            cout << "r = " << r << endl;
         getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
//         if (index[0] == 24 && index[1] == 28 && index[2] == 39)
//            cout << "alpha = " << alpha << endl;
         beta = 1.0-alpha;
         ebar = alpha*ethere+beta*ehere;
         for (j = 0; j < grid.dim; j++)
            M[r][j] = -gamma[r][(s+1)/2]*gamma[j][(s+1)/2]*beta*(ehere-ethere)/ebar*
                       tangent[r]*tangent[j];
         M[r][r] += 1.0;
            
         for (j = 0; j < grid.dim; j++)
            f[r][j] = (gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                       (1-gamma[j][(s+1)/2])*s*tangent[j])/(grid.dx[r]*ebar);
         f[r][r] += (s*ethere)/(grid.dx[r]*ebar);
         c[r] = s*ethere;
         for (j = 0; j < grid.dim; j++)
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
      for (r = 0; r < grid.dim; r++)
      {
         for (j = 0; j < grid.dim; j++)
            temp[j] = f[j][r];
         forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
         rindex[r] = index[r]+s;
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
   free_matrix(f,grid.dim-1,grid.dim-1);
}

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




void perturb(double ***S, double tol, GridData &grid)
{
   int tindex[grid.dim];
   int i;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      perturbelt(S,tindex,tol);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}

void perturb(double ***S, double tol, PBData &pb, GridData &grid)
{
   int tindex[grid.dim];
   int i;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      perturbelt(S,tindex,tol,pb);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}
// pertube for maxstep number of steps according to status
void perturbstatus(double ***S, double tol, int maxsteps, PBData &pb, GridData &grid)
{
   int i, step, numperturb, tindex[grid.dim];
   char yesperturb = 1;

   for (step = 1; step <= maxsteps && yesperturb; step++)
   {
      numperturb = 0;
      yesperturb = 0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         if (perturbstatuselt(S,tindex,tol,pb,grid))
         {
            yesperturb = 1;
            numperturb++;
         }
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      cout << "Perturbed " << numperturb << " points at step " << step << endl;
   }
}
// if surface on grid point, within tol, shift to tol
void perturbelt(double ***S, int *index, double tol)
{
   if (evalarray(S,index) >= 0.0 && evalarray(S,index) < tol)
      setvalarray(S,index,tol);
   else if (evalarray(S,index) < 0.0 && evalarray(S,index) > -tol)
      setvalarray(S,index,-tol);
}
// if outside, within tol, push to tol depending on epsilon
void perturbelt(double ***S, int *index, double tol, PBData &pb)
{
   if (evalarray(S,index) >= 0.0 && evalarray(S,index) < tol)
   {
      if (pb.epsilonp >= pb.epsilonm)
         setvalarray(S,index,tol);
      else
         setvalarray(S,index,-tol);
   }
   else if (evalarray(S,index) < 0.0 && evalarray(S,index) > -tol)
   {
      if (pb.epsilonm >= pb.epsilonp)
         setvalarray(S,index,-tol);
      else
         setvalarray(S,index,tol);
   }
}
// perturb according to status
char perturbstatuselt(double ***S, int *index, double tol, PBData &pb, GridData &grid)
{

   int r, s, t, thesign, count[2], rindex[grid.dim];

   if (fabs(evalarray(S,index)) <= tol)// first perturbe according to epsilon, move point to larger epsilon side
   {  
      if (pb.epsilonp >= pb.epsilonm){
       setvalarray(S,index,tol);
      }
      else{
       setvalarray(S,index,-tol);
      }
   }
   // if index is cim3 points,  count the number of cim points in nbr before and after flipping
   // if flipping reduce cim points, flip back (prefer more cim points) and return 1
   if (checkcim3(S,index,grid) == 3) 
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];
      for (t = 0; t <= 1; t++)
      {
         count[t] = 0;
         for (r = 0; r < grid.dim; r++)
         {
            for (s = -1; s <= 1; s += 2)
            {
               rindex[r] = index[r]+s;
               if (checkcim3(S,rindex,grid) == 3)
                  (count[t])++;
            }
            rindex[r] = index[r];
         }
         setvalarray(S,index,-evalarray(S,index));
      }
      if (count[1] > count[0]) // original count[1] < count[0]
      {
         setvalarray(S,index,-evalarray(S,index));
         printf("Perturb (%d %d %d) to %.6e\n", index[0],index[1],index[2],evalarray(S,index));
         return 1;
      }
   }

   return 0;
}



void Atransposesmall(SparseElt2**** &B, SparseElt2**** &A, GridData &grid, double ***S,
                     PBData &pb)
{
   int i, r, m, n, tindex[grid.dim], sindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   double ehere, value;

// copy A to B and order B
   clearsparse(B,grid);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (evalarray(A,tindex) == NULL)
      {
// create rows in B when they don't exist in A
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;

         value = 0.0;
         for (m = 0; m < grid.dim; m++)
            value += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
         sparseorder(tindex,tindex,B,value,grid);

         for (r = 0; r < grid.dim; r++) 
            sindex[r] = tindex[r];
         for (m = 0; m < grid.dim; m++)
         {
            for (n = -1; n <= 1; n += 2)
            {
               sindex[m] = tindex[m]+n;
               if (sindex[m] >= 0 && sindex[m] <= grid.nx[m])
               {
                  value = -ehere/(grid.dx[m]*grid.dx[m]);
                  sparseorder(sindex,tindex,B,value,grid);
               }
            }
            sindex[m] = tindex[m];
         }
      }
      else
         for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
            sparseorder((*current).cindex,tindex,B,(*current).val,grid);
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}

void copyfromsmall(SparseElt2**** &B, SparseElt2**** &A, GridData &grid, double ***S, 
                   PBData &pb)
{
   int i, r, m, n, tindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   double ehere;

// copy A to B and order B
   clearsparse(B,grid);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(B,tindex,NULL);
      if (evalarray(A,tindex) == NULL)
      {
// create rows in B when they don't exist in A
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
         temp = new SparseElt2;
         (*temp).cindex = new int[grid.dim];
         for (r = 0; r < grid.dim; r++) 
            (*temp).cindex[r] = tindex[r];
         (*temp).val = 0.0;
         for (m = 0; m < grid.dim; m++)
            (*temp).val += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
         (*temp).next = NULL;
         setvalarray(B,tindex,temp);

         for (m = 0; m < grid.dim; m++)
            for (n = -1; n <= 1; n += 2)
            {
               temp = new SparseElt2;
               (*temp).cindex = new int[grid.dim];
               for (r = 0; r < grid.dim; r++) 
                  (*temp).cindex[r] = tindex[r];
               (*temp).cindex[m] = tindex[m]+n;
               if ((*temp).cindex[m] >= 0 && (*temp).cindex[m] <= grid.nx[m])
               {
                  (*temp).val = -ehere/(grid.dx[m]*grid.dx[m]);
                  for (prev = NULL,current2 = evalarray(B,tindex); 
                       current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                           sub2ind((*temp).cindex,grid.nx,grid.dim); 
                       prev = current2,current2 = (*current2).next);
                  (*temp).next = current2;
                  if (prev != NULL)
                     (*prev).next = temp;
                  else
                     setvalarray(B,tindex,temp);
               }
            }
      }
      else
         for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
         {
            temp = new SparseElt2;
            (*temp).cindex = new int[grid.dim];
            for (r = 0; r < grid.dim; r++) 
               (*temp).cindex[r] = (*current).cindex[r];
            (*temp).val = (*current).val;
            for (prev = NULL,current2 = evalarray(B,tindex); 
                 current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                     sub2ind((*temp).cindex,grid.nx,grid.dim); 
                 prev = current2,current2 = (*current2).next);
            (*temp).next = current2;
            if (prev != NULL)
               (*prev).next = temp;
            else
               setvalarray(B,tindex,temp);
         }   
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}

void copyfromsmall(SparseElt2**** &B, SparseElt2**** &A, GridData &grid, double ***a, 
                   double ***S, PBData &pb)
{
   int i, r, m, n, tindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   double ehere;

// copy A to B and order B
   clearsparse(B,grid);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(B,tindex,NULL);
      if (evalarray(A,tindex) == NULL)
      {
// create rows in B when they don't exist in A
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
         temp = new SparseElt2;
         (*temp).cindex = new int[grid.dim];
         for (r = 0; r < grid.dim; r++) 
            (*temp).cindex[r] = tindex[r];
         (*temp).val = evalarray(a,tindex);
         for (m = 0; m < grid.dim; m++)
            (*temp).val += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
         (*temp).next = NULL;
         setvalarray(B,tindex,temp);

         for (m = 0; m < grid.dim; m++)
            for (n = -1; n <= 1; n += 2)
            {
               temp = new SparseElt2;
               (*temp).cindex = new int[grid.dim];
               for (r = 0; r < grid.dim; r++) 
                  (*temp).cindex[r] = tindex[r];
               (*temp).cindex[m] = tindex[m]+n;
               if ((*temp).cindex[m] >= 0 && (*temp).cindex[m] <= grid.nx[m])
               {
                  (*temp).val = -ehere/(grid.dx[m]*grid.dx[m]);
                  for (prev = NULL,current2 = evalarray(B,tindex); 
                       current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                           sub2ind((*temp).cindex,grid.nx,grid.dim); 
                       prev = current2,current2 = (*current2).next);
                  (*temp).next = current2;
                  if (prev != NULL)
                     (*prev).next = temp;
                  else
                     setvalarray(B,tindex,temp);
               }
            }
      }
      else
         for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
         {
            temp = new SparseElt2;
            (*temp).cindex = new int[grid.dim];
            for (r = 0; r < grid.dim; r++) 
               (*temp).cindex[r] = (*current).cindex[r];
            (*temp).val = (*current).val;
            for (prev = NULL,current2 = evalarray(B,tindex); 
                 current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                     sub2ind((*temp).cindex,grid.nx,grid.dim); 
                 prev = current2,current2 = (*current2).next);
            (*temp).next = current2;
            if (prev != NULL)
               (*prev).next = temp;
            else
               setvalarray(B,tindex,temp);
         }   
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}

void ILUsmall(SparseElt2**** &M, SparseElt2**** &A, GridData &grid, double ***S, 
              PBData &pb)
{
   int i, r, m, n, tindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   double ehere;
   int N = 1;
   for (r = 0; r < grid.dim; r++)
      N *= grid.nx[r]+1;

   copyfromsmall(M,A,grid,S,pb);
/*
// copy A to M and order M
   clearsparse(M,grid);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(M,tindex,NULL);
      if (evalarray(A,tindex) == NULL)
      {
// create rows in M when they don't exist in A
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
         temp = new SparseElt2;
         (*temp).cindex = new int[grid.dim];
         for (r = 0; r < grid.dim; r++) 
            (*temp).cindex[r] = tindex[r];
         (*temp).val = 0.0;
         for (m = 0; m < grid.dim; m++)
            (*temp).val += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
         (*temp).next = NULL;
         setvalarray(M,tindex,temp);

         for (m = 0; m < grid.dim; m++)
            for (n = -1; n <= 1; n += 2)
            {
               temp = new SparseElt2;
               (*temp).cindex = new int[grid.dim];
               for (r = 0; r < grid.dim; r++) 
                  (*temp).cindex[r] = tindex[r];
               (*temp).cindex[m] = tindex[m]+n;
               if ((*temp).cindex[m] >= 0 && (*temp).cindex[m] <= grid.nx[m])
                  (*temp).val = -ehere/(grid.dx[m]*grid.dx[m]);
               for (prev = NULL,current2 = evalarray(M,tindex); 
                    current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                        sub2ind((*temp).cindex,grid.nx,grid.dim); 
                    prev = current2,current2 = (*current2).next);
               (*temp).next = current2;
               if (prev != NULL)
                  (*prev).next = temp;
               else
                  setvalarray(M,tindex,temp);
            }
      }
      else
         for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
         {
            temp = new SparseElt2;
            (*temp).cindex = new int[grid.dim];
            for (r = 0; r < grid.dim; r++) 
               (*temp).cindex[r] = (*current).cindex[r];
            (*temp).val = (*current).val;
            for (prev = NULL,current2 = evalarray(M,tindex); 
                 current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                     sub2ind((*temp).cindex,grid.nx,grid.dim); 
                 prev = current2,current2 = (*current2).next);
            (*temp).next = current2;
            if (prev != NULL)
               (*prev).next = temp;
            else
               setvalarray(M,tindex,temp);
         }   
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
*/

// work on M
   for (i = 0; i < N; i++)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
// row is tindex, column is (*current).cindex
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
      {
// (*current2).cindex is different columns of row tindex
         for (current2 = evalarray(M,tindex); current2 != current && 
              sub2ind((*current2).cindex,grid.nx,grid.dim) < 
              min(sub2ind(tindex,grid.nx,grid.dim),
                  sub2ind((*current).cindex,grid.nx,grid.dim)); 
              current2 = (*current2).next)
         {
// (*temp).cindex is different columns of row (*current2).cindex
            for (temp = evalarray(M,(*current2).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val -= ((*current2).val)*((*temp).val);
         }
// if working on L
         if (sub2ind((*current).cindex,grid.nx,grid.dim) < 
             sub2ind(tindex,grid.nx,grid.dim))
         {
            for (temp = evalarray(M,(*current).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val /= (*temp).val;
         }
      }
   }
}

void ILUsmall(SparseElt2**** &M, SparseElt2**** &A, GridData &grid, double ***a, 
              double ***S, PBData &pb)
{
   int i, r, m, n, tindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   double ehere;
   int N = 1;
   for (r = 0; r < grid.dim; r++)
      N *= grid.nx[r]+1;

   copyfromsmall(M,A,grid,a,S,pb);

// work on M
   for (i = 0; i < N; i++)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
// row is tindex, column is (*current).cindex
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
      {
// (*current2).cindex is different columns of row tindex
         for (current2 = evalarray(M,tindex); current2 != current && 
              sub2ind((*current2).cindex,grid.nx,grid.dim) < 
              min(sub2ind(tindex,grid.nx,grid.dim),
                  sub2ind((*current).cindex,grid.nx,grid.dim)); 
              current2 = (*current2).next)
         {
// (*temp).cindex is different columns of row (*current2).cindex
            for (temp = evalarray(M,(*current2).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val -= ((*current2).val)*((*temp).val);
         }
// if working on L
         if (sub2ind((*current).cindex,grid.nx,grid.dim) < 
             sub2ind(tindex,grid.nx,grid.dim))
         {
            for (temp = evalarray(M,(*current).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val /= (*temp).val;
         }
      }
   }
}

void leftmultILUinv(double ***y, SparseElt2**** &M, double ***x, GridData &grid)
{
   int i, tindex[grid.dim];
   double d;
   SparseElt2 *current;
   int N = 1;
   for (i = 0; i < grid.dim; i++)
      N *= grid.nx[i]+1;
/*
   clock_t cstart1, cend1;
   clock_t cstart2, cend2;
   clock_t cstart3, cend3;
   clock_t cstart4, cend4;
   clock_t cstart5, cend5;
   clock_t cstart6, cend6;
   clock_t cstart7, cend7;
   clock_t cstart8, cend8;
   double temp;
   
   cstart1 = clock ();
   for (i = 0; i < N; i++)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
      setvalarray(y,tindex,0.0);
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
         setvalarray(y,tindex,evalarray(y,tindex)+
                              (*current).val*evalarray(x,(*current).cindex));
   }
   cend1 = clock ();
   cout << "mx-vec mult time = " << (double) (cend1-cstart1)/CLOCKS_PER_SEC << endl;
   cstart3 = clock ();
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(y,tindex,0.0);
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
         setvalarray(y,tindex,evalarray(y,tindex)+
                              (*current).val*evalarray(x,(*current).cindex));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend3 = clock ();
   cout << "mx-vec mult time = " << (double) (cend3-cstart3)/CLOCKS_PER_SEC << endl;
   cstart4 = clock ();
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = 0.0;
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
         temp += (*current).val*evalarray(x,(*current).cindex);
      setvalarray(y,tindex,temp);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend4 = clock ();
   cout << "mx-vec mult time = " << (double) (cend4-cstart4)/CLOCKS_PER_SEC << endl;
   cstart6 = clock ();
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(x,tindex);
      for (current = evalarray(M,tindex); 
           current != NULL && 
           ((*current).cindex[0] < tindex[0] || 
            ((*current).cindex[0] == tindex[0] && (*current).cindex[1] < tindex[1]) || 
            ((*current).cindex[0] == tindex[0] && (*current).cindex[1] == tindex[1] &&
             (*current).cindex[2] < tindex[2]));
           current = (*current).next)
         temp -= (*current).val*evalarray(y,(*current).cindex);
      setvalarray(y,tindex,temp);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend6 = clock ();
   cout << "forward solve time = " << (double) (cend6-cstart6)/CLOCKS_PER_SEC << endl;
   cstart7 = clock ();
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(x,tindex);
      for (current = evalarray(M,tindex); 
           current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) < 
                              sub2ind(tindex,grid.nx,grid.dim);
           current = (*current).next)
         temp -= (*current).val*evalarray(y,(*current).cindex);
      setvalarray(y,tindex,temp);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend7 = clock ();
   cout << "forward solve time = " << (double) (cend7-cstart7)/CLOCKS_PER_SEC << endl;
   int r;
   cstart8 = clock ();
   r = 0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(x,tindex);
      for (current = evalarray(M,tindex); 
           current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) < r;
           current = (*current).next)
         temp -= (*current).val*evalarray(y,(*current).cindex);
      setvalarray(y,tindex,temp);
      r++;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend8 = clock ();
   cout << "forward solve time = " << (double) (cend8-cstart8)/CLOCKS_PER_SEC << endl;
*/
//   cstart2 = clock ();
   for (i = 0; i < N; i++)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
      setvalarray(y,tindex,evalarray(x,tindex));
      for (current = evalarray(M,tindex); 
           current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) < i;
           current = (*current).next)
         setvalarray(y,tindex,evalarray(y,tindex)-
                              (*current).val*evalarray(y,(*current).cindex));
   }
//   cend2 = clock ();
//   cout << "forward solve time = " << (double) (cend2-cstart2)/CLOCKS_PER_SEC << endl;
//   cstart5 = clock ();
   for (i = N-1; i >= 0; i--)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
      for (current = evalarray(M,tindex); 
           current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) < i;
           current = (*current).next);
      d = (*current).val;
      for (current = (*current).next; current != NULL; current = (*current).next)
         setvalarray(y,tindex,evalarray(y,tindex)-
                              (*current).val*evalarray(y,(*current).cindex));
      setvalarray(y,tindex,evalarray(y,tindex)/d);
   }
//   cend5 = clock ();
//   cout << "backward solve time = " << (double) (cend5-cstart5)/CLOCKS_PER_SEC << endl;
//   exit(1);
}

double getpsivn(double ***u, double ***S, int *index, int rstar, PBData &pb, 
                GridData grid)
{
// gets normal velocity at interface point neighboring location index in the 
// positive direction of the rstar dimension
   int r;
   int rindex[grid.dim];
   double grad[grid.dim], normal[grid.dim], tangent[grid.dim];
   double alpha, ehere, ethere, thesign, v1, v2, vn, temp;

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

   double v1e, v2e;
   int stat1, stat2;
   double uint, uint2;
   int tindex[grid.dim];
   for (r = 0; r < grid.dim; r++)
      tindex[r] = index[r];
   tindex[rstar] = index[rstar]-1;

   uint = getinterfacegrad3(grad,u,S,index,rstar,1,pb,grid);
   cout << uint << endl;
   cout << evalarray(u,index) << " " << evalarray(S,index) << endl;
   cout << evalarray(u,tindex) << " " << evalarray(S,tindex) << endl;
   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   rindex[rstar] = index[rstar]+1;
   getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
   temp = 0.0;
   for (r = 0; r < grid.dim; r++)
      temp += grad[r]*normal[r];
   v1 = thesign*ehere*temp*temp;
   v1e = thesign*ehere*ethere*ethere;
   temp = 0.0;
   for (r = 0; r < grid.dim; r++)
      temp += grad[r]*grad[r];
   v2 = 0.5*thesign*ehere*temp;
   v2e = 0.5*thesign*ehere*4.0*ethere*ethere*0.25;
   stat1 = getstatus3(S,index,grid);
   cout << stat1 << endl;

   for (r = 0; r < grid.dim; r++)
      tindex[r] = rindex[r];
   tindex[rstar] = rindex[rstar]+1;
   uint2 = getinterfacegrad3(grad,u,S,rindex,rstar,-1,pb,grid);
   uint = (uint+uint2)/2.0;
   cout << uint2 << endl;
   cout << evalarray(u,rindex) << " " << evalarray(S,rindex) << endl;
   cout << evalarray(u,tindex) << " " << evalarray(S,tindex) << endl;
   temp = 0.0;
   for (r = 0; r < grid.dim; r++)
      temp += grad[r]*normal[r];
   v1 -= thesign*ethere*temp*temp;
   v1e -= thesign*ethere*ehere*ehere;
   temp = 0.0;
   for (r = 0; r < grid.dim; r++)
      temp += grad[r]*grad[r];
   v2 -= 0.5*thesign*ethere*temp;
   v2e -= 0.5*thesign*ethere*4.0*ehere*ehere*0.25;
   stat2 = getstatus3(S,rindex,grid);
   cout << stat2 << endl;

   vn = v1-v2;
   vn += -Bval(uint,pb);
//   getchar();

   return vn;
}

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

   uint = getinterfacegrad3(grad,u,S,tindex,rstar,rindex[rstar]-tindex[rstar],pb,grid);
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
   return getpsivn(psi,S,index,rstar,pb,grid)+getLJvn(LJ,S,index,rstar,pb,grid);
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

   getinterfaceDu(uhere,Duhere,index,rstar,sstar,u,S,pb,grid);
   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   rindex[rstar] = index[rstar]+sstar;
   getinterfaceDu(uthere,Duthere,rindex,rstar,-sstar,u,S,pb,grid);

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

//   double maxuinterr = 0.0, uhere2, uthere2, Duhere2[grid.dim], Duthere2[grid.dim], 
//          maxDuerr[grid.dim];

//   for (r = 0; r < grid.dim; r++)
//      maxDuerr[r] = 0.0;
//   getinterfaceDu(uhere2,Duhere2,index,rstar,sstar,u,S,pb,grid);
//   evalfromstorage(uhere,Duhere,index,rstar,sstar,mid,Dusmall,smallsize,u,S,pb,grid);
   evalfromstorage(uhere,Duhere,index,rstar,sstar,Dusmall,smallsize,u,S,pb,grid);
//   if (fabs(uhere2-uhere) > maxuinterr)
//      maxuinterr = fabs(uhere2-uhere);
//   for (r = 0; r < grid.dim; r++)
//      if (fabs(Duhere2[r]-Duhere[r]) > maxDuerr[r])
//         maxDuerr[r] = fabs(Duhere2[r]-Duhere[r]);
   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   rindex[rstar] = index[rstar]+sstar;
//   getinterfaceDu(uthere2,Duthere2,rindex,rstar,-sstar,u,S,pb,grid);
//   evalfromstorage(uthere,Duthere,rindex,rstar,-sstar,mid,Dusmall,smallsize,u,S,pb,
//                   grid);
   evalfromstorage(uthere,Duthere,rindex,rstar,-sstar,Dusmall,smallsize,u,S,pb,grid);
//   if (fabs(uthere2-uthere) > maxuinterr)
//      maxuinterr = fabs(uthere2-uthere);
//   for (r = 0; r < grid.dim; r++)
//      if (fabs(Duthere2[r]-Duthere[r]) > maxDuerr[r])
//         maxDuerr[r] = fabs(Duthere2[r]-Duthere[r]);
//   if (maxuinterr > 1.0e-11 || maxDuerr[0] > 1.0e-11 || maxDuerr[1] > 1.0e-11 || 
//       maxDuerr[2] > 1.0e-11)
//   {
//      cout << "maxuinterr = " << maxuinterr << endl;
//      cout << "maxDuerr = " << maxDuerr[0] << " " << maxDuerr[1] << " " << maxDuerr[2] 
//           << endl;
//      getchar();
//   }

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
/*
      cout << grid.maxvn << endl;
      cout << fabs(Duhere[0]-getDu(index,0,rstar,sstar,alpha,thesign,grid)) << " "
           << fabs(Duhere[1]-getDu(index,1,rstar,sstar,alpha,thesign,grid)) << " "
           << fabs(Duhere[2]-getDu(index,2,rstar,sstar,alpha,thesign,grid)) << endl;
      cout << fabs(Duthere[0]-getDu(index,0,rstar,sstar,alpha,-thesign,grid)) << " "
           << fabs(Duthere[1]-getDu(index,1,rstar,sstar,alpha,-thesign,grid)) << " "
           << fabs(Duthere[2]-getDu(index,2,rstar,sstar,alpha,-thesign,grid)) << endl;
      cout << vn << " " << (getDu(index,0,rstar,sstar,alpha,1.0,grid)-
                            getDu(index,0,rstar,sstar,alpha,-1.0,grid))*normal[0]+
                           (getDu(index,1,rstar,sstar,alpha,1.0,grid)-
                            getDu(index,1,rstar,sstar,alpha,-1.0,grid))*normal[1]+
                           (getDu(index,2,rstar,sstar,alpha,1.0,grid)-
                            getDu(index,2,rstar,sstar,alpha,-1.0,grid))*normal[2] << " "
           << " " << 2.0*(1.0-pb.epsilonp/pb.epsilonm)*0.5*
                     exp(2.0*(1.0-pb.epsilonp/pb.epsilonm)*grid.t) << endl;
      getchar();
*/
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

void outputstatus(ofstream& outfile, double ***phi, GridData &grid)
{
   int i, tindex[grid.dim];

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (evalarray(phi,tindex) < 0.0)
         outfile << -getstatus3(phi,tindex,grid) << " ";
      else
         outfile << getstatus3(phi,tindex,grid) << " ";
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
         if (i == grid.dim-1)
            outfile << endl;
      }
   }
}

void cim2again3(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
                PBData &pb, GridData &grid)
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
   double sigma, d[grid.dim], thesign, tau, Dtau[grid.dim];

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

   if (index[0] == 25 && index[1] == 25 && index[2] == 16)
      cout << "HERE" << endl;
   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      tindex[r] = 2;
   for (r = 0; r < grid.dim; r++)
   {
      rindex[r] = index[r]+sk[r];
      getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
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
            tindex[r] = 2+(s-1)*sk[r];
            setvalarray(f[r],tindex,evalarray(f[r],tindex)+a[s]);
         }
         tindex[r] = 2-sk[r];
         setvalarray(f[r],tindex,evalarray(f[r],tindex)-sk[r]*C*tangent[r]*sk[r]);
         tindex[r] = 2;
         setvalarray(f[r],tindex,evalarray(f[r],tindex)+sk[r]*C*tangent[r]*sk[r]);

         M[r][r] = 1.0-abs(sk[r])*(0.5+alpha)*bk*tangent[r]*tangent[r];

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

               if (abs(tmps) == 1)
               {
                  tindex[j] = 2-tmps;
                  rindex[j] = index[j]-tmps;
                  if (evalarray(S,rindex)*Sval <= 0.0)
                     cout << "wrong 4" << endl;
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                          sk[r]*C*tangent[j]*(1.0+alpha)*tmps);
                  tindex[r] = 2-sk[r];
                  rindex[r] = index[r]-sk[r];
                  if (evalarray(S,rindex)*Sval <= 0.0)
                     cout << "wrong 5" << endl;
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                          sk[r]*C*tangent[j]*alpha*tmps);
                  tindex[j] = 2;
                  rindex[j] = index[j];
                  if (evalarray(S,rindex)*Sval <= 0.0)
                     cout << "wrong 6" << endl;
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                          sk[r]*C*tangent[j]*alpha*tmps);
                  tindex[r] = 2;
                  rindex[r] = index[r];
                  if (evalarray(S,rindex)*Sval <= 0.0)
                     cout << "wrong 7" << endl;
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                          sk[r]*C*tangent[j]*(1.0+alpha)*tmps);
               }
               else
               {
                  for (s = -1; s <= 1; s += 2)
                  {
                     tindex[j] = 2+s;
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)
                        cout << "wrong 8" << endl;
                     setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                             sk[r]*C*tangent[j]*(1.0+alpha)*0.5*s);
                  }
                  tindex[r] = 2-sk[r];
                  rindex[r] = index[r]-sk[r];
                  for (s = -1; s <= 1; s += 2)
                  {
                     tindex[j] = 2+s;
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)
                        cout << "wrong 9" << endl;
                     setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                             sk[r]*C*tangent[j]*alpha*0.5*s);
                  }
                  tindex[r] = 2;
                  tindex[j] = 2;
                  rindex[r] = index[r];
                  rindex[j] = index[j];
               }
            }

// ADDING
         getsigma(sigma,index,r,sk[r],alpha,normal,pb,grid);
         gettau(tau,index,r,sk[r],alpha,grid);
         getDtau(Dtau,index,r,sk[r],alpha,grid);
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
            tindex[r] = 2+s;
            setvalarray(f[r],tindex,evalarray(f[r],tindex)+1.0);
         }
         tindex[r] = 2;
         setvalarray(f[r],tindex,evalarray(f[r],tindex)-2.0);
      }
   }

   gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);
    
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] < 5)
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
            rindex[r] = index[r]-2+tindex[r];
         sparse2(index,rindex,A,value,grid);
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] >= 5; i--)
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

   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);
   free_matrix(f,grid.dim-1,4,4,4);
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

void cim3(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
          PBData &pb, GridData &grid)
{
   int r, s, t, i, m, n;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];
   double d0[grid.dim], b0[2*grid.dim], c0[2*grid.dim]; 
   double ***dcoef[grid.dim], ***bcoef[2*grid.dim];
//   double ***D1[grid.dim], ***D2[grid.dim][grid.dim];
   double ***D1[grid.dim], ***D2[grid.dim][3];
   double **LU, **G, **B, **Dn, value, Sval; 
   double Dndott[grid.dim], Dndots[grid.dim], dotDndot[grid.dim];
   int sk2[grid.dim][grid.dim][4];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double alpha, beta, tangent1[grid.dim], tangent2[grid.dim], normal[grid.dim];
   double ethere, ehere, jumpfe;
   double temp[2*grid.dim];
   int sk[grid.dim];
   double sigma, tau, Dsigma[grid.dim], thesign;
   double Dtau[grid.dim], **D2tau, dotD2taudot[grid.dim], Dtaudot[grid.dim];

   double realb[2*grid.dim], realjumpuxx[2*grid.dim], realux[grid.dim], 
          dotrealuxxdot[2], realvalue, approxux[grid.dim], reald[grid.dim];

   cout << "in cim3" << endl;
   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   G = matrix(grid.dim-1,grid.dim-1);
   B = matrix(2*grid.dim-1,2*grid.dim-1);
   Dn = matrix(grid.dim-1,grid.dim-1);
   D2tau = matrix(grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      D1[r] = matrix(2,2,2);
      for (s = 0; s < grid.dim; s++)
         D2[r][s] = matrix(2,2,2);
      dcoef[r] = matrix(2,2,2);
   }
   for (r = 0; r < 2*grid.dim; r++)
      bcoef[r] = matrix(2,2,2);

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

   cout << "index = " << index[0] << " " << index[1] << " " << index[2] << endl;
   setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid));
   cout << "here" << endl;

   for (r = 0; r < grid.dim; r++)
      sk[r] = gamma[r][1]-gamma[r][0];
   cout << "here" << endl;

   int j;
   double jumpuxx[grid.dim], uxy[grid.dim][grid.dim], ux[grid.dim];
   double x[grid.dim];
   double ***ucoef[grid.dim], **uxxcoef;
   uxxcoef = matrix(grid.dim,grid.dim);
   for (r = 0; r < grid.dim; r++)
      ucoef[r] = matrix(2,2,2);
   char theorder = 2;
   double jumpu, ***jumpucoef, jumpuxxcoef[grid.dim];
   jumpucoef = matrix(2,2,2);

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
   {
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 0;
      while (sindex[0] < 3)
      {
         setvalarray(dcoef[r],sindex,0.0);
         for (s = 0; s < 2*grid.dim; s++)
            setvalarray(bcoef[s],sindex,0.0);
         for (s = 0; s < grid.dim; s++)
         {
            setvalarray(D1[s],sindex,0.0);
            for (n = 0; n < grid.dim; n++)
               setvalarray(D2[s][n],sindex,0.0);
         }
         (sindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
         {
            sindex[i] = 0;
            (sindex[i-1])++;
         }
      }
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 1;

      cout << "r = " << r << ", sk[r] = " << sk[r] << endl;
      if (abs(sk[r]) == 1)
      {
         rindex[r] = index[r]+sk[r];
         for (m = 0; m < grid.dim; m++)
         {
            if (sk[m] == 0)
            {
               for (s = -1; s <= 1; s += 2)
               {
                  sindex[m] = 1+s;
                  setvalarray(D1[m],sindex,s/(2.0*grid.dx[m]));
               }
               sindex[m] = 1;
            }
            else
            {
               for (s = 0; s <= 1; s++)
               {
                  sindex[m] = 1-s*sk[m];
                  setvalarray(D1[m],sindex,(1-2*s)*sk[m]/grid.dx[m]);
               }
               sindex[m] = 1;
            }
         }
         cout << "created D1" << endl;
         for (m = 0; m < grid.dim; m++)
            for (n = m+1; n < grid.dim; n++)
            {
// get sk2
               getsk2(sk2[m][n],m,n,index,S,grid);
               for (s = 0; s < 2; s++)
               {
                  sindex[m] = 1+sk2[m][n][s];
                  for (t = 0; t < 2; t++)
                  {
                     sindex[n] = 1+sk2[m][n][2+t];
                     setvalarray(D2[m][n],sindex,(2*s-1)*(2*t-1)/
                                                 ((sk2[m][n][1]-sk2[m][n][0])*
                                                  (sk2[m][n][3]-sk2[m][n][2])*
                                                  grid.dx[m]*grid.dx[n]));
                  }
                  sindex[n] = 1;
               }
               sindex[m] = 1;
            }
         for (j = 0; j < grid.dim; j++)
            tindex[j] = index[j];
         for (m = 0; m < grid.dim; m++)
         {
            ux[m] = 0.0;
            for (n = m+1; n < grid.dim; n++)
               uxy[m][n] = 0.0;
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            for (j = 0; j < grid.dim; j++)
               tindex[j] = index[j]-1+sindex[j];
            for (m = 0; m < grid.dim; m++)
            {
               ux[m] += evalarray(D1[m],sindex)*getu(tindex,0,0,0.0,thesign,grid);
               for (n = m+1; n < grid.dim; n++)
                  uxy[m][n] += evalarray(D2[m][n],sindex)*getu(tindex,0,0,0.0,thesign,
                                                               grid);
            }
            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (j = 0; j < grid.dim; j++)
            tindex[j] = index[j];
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;
         for (m = 0; m < grid.dim; m++)
            cout << m << " " << ux[m] << " " << getDu(index,m,0,0,0.0,thesign,grid)
                 << endl;
         for (m = 0; m < grid.dim; m++)
            for (n = m+1; n < grid.dim; n++)
               cout << m << " " << n << " " << uxy[m][n] << " " 
                    << getD2u(index,m,n,0,0,0.0,thesign,grid) << " " 
                    << sk2[m][n][0] << " " << sk2[m][n][1] << " " << sk2[m][n][2] 
                    << " " << sk2[m][n][3] << endl;

// get Dn, Dsigma, jumpfe
// problem if tangent1 = 0
         getinterfaceinfo(alpha,tangent1,tangent2,normal,sigma,Dn,Dsigma,jumpfe,S,index,
                          rindex,pb,grid);
         beta = 1.0-alpha;
//         cout << index[0] << " " << index[1] << " " << index[2] << endl;
//         cout << rindex[0] << " " << rindex[1] << " " << rindex[2] << endl;
//         cout << evalarray(S,index) << " " << evalarray(S,rindex) << endl;
//         cout << "alpha = " << alpha << endl;
//         cout << "beta = " << beta << endl;
         cout << tangent1[0] << " " << tangent1[1] << " " << tangent1[2] << endl;
         cout << tangent2[0] << " " << tangent2[1] << " " << tangent2[2] << endl;
         cout << normal[0] << " " << normal[1] << " " << normal[2] << endl;
//         cout << "sigma = " << sigma << endl;
//         cout << "Dsigma = " << Dsigma[0] << " " << Dsigma[1] << " " << Dsigma[2] 
//              << endl;
//         cout << "jumpfe = " << jumpfe << endl;
//         for (m = 0; m < grid.dim; m++)
//         {
//            for (n = 0; n < grid.dim; n++)
//               cout << Dn[m][n] << " ";
//            cout << endl;
//         }
         cout << "check " << -pb.epsilonp*(getD2u(index,0,0,0,0,0.0,1,grid)+
                                           getD2u(index,1,1,0,0,0.0,1,grid)+
                                           getD2u(index,2,2,0,0,0.0,1,grid))-
                             getf(index,0,0,0.0,1,pb,grid) << " "
                          << -pb.epsilonm*(getD2u(index,0,0,0,0,0.0,-1,grid)+
                                           getD2u(index,1,1,0,0,0.0,-1,grid)+
                                           getD2u(index,2,2,0,0,0.0,-1,grid))-
                             getf(index,0,0,0.0,-1,pb,grid) << endl;
         cout << "check " << pb.epsilonp*(getDu(index,0,r,sk[r],alpha,1,grid)*normal[0]+
                                          getDu(index,1,r,sk[r],alpha,1,grid)*normal[1]+
                                          getDu(index,2,r,sk[r],alpha,1,grid)*normal[2])-
                             pb.epsilonm*(getDu(index,0,r,sk[r],alpha,-1,grid)*normal[0]+
                                          getDu(index,1,r,sk[r],alpha,-1,grid)*normal[1]+
                                          getDu(index,2,r,sk[r],alpha,-1,grid)*normal[2])-
                             sigma << endl;
         for (j = 0; j < grid.dim; j++)
            cout << "check " << pb.epsilonp*(getD2u(index,j,0,r,sk[r],alpha,1,grid)*
                                             normal[0]+
                                             getD2u(index,j,1,r,sk[r],alpha,1,grid)*
                                             normal[1]+
                                             getD2u(index,j,2,r,sk[r],alpha,1,grid)*
                                             normal[2]+
                                             getDu(index,0,r,sk[r],alpha,1,grid)*
                                             Dn[0][j]+
                                             getDu(index,1,r,sk[r],alpha,1,grid)*
                                             Dn[1][j]+
                                             getDu(index,2,r,sk[r],alpha,1,grid)*
                                             Dn[2][j])-
                                pb.epsilonm*(getD2u(index,j,0,r,sk[r],alpha,-1,grid)*
                                             normal[0]+
                                             getD2u(index,j,1,r,sk[r],alpha,-1,grid)*
                                             normal[1]+
                                             getD2u(index,j,2,r,sk[r],alpha,-1,grid)*
                                             normal[2]+
                                             getDu(index,0,r,sk[r],alpha,-1,grid)*
                                             Dn[0][j]+
                                             getDu(index,1,r,sk[r],alpha,-1,grid)*
                                             Dn[1][j]+
                                             getDu(index,2,r,sk[r],alpha,-1,grid)*
                                             Dn[2][j])-
                                Dsigma[j] << endl;
         sub2coord(x,index,grid);
         x[r] += sk[r]*alpha*grid.dx[r];
         for (m = 0; m < grid.dim; m++)
            for (n = 0; n < grid.dim; n++)
            {
               cout << "check " << m << " " << n << " " << Dn[m][n];
               if (m != n)
                  cout << " " << -x[m]*x[n]/((x[0]*x[0]+x[1]*x[1]+x[2]*x[2])*
                                             sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])) << endl;
               else
                  cout << " " << (1.0-x[m]*x[n]/(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/
                                 sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << endl;
            }
//         cout << "   compare getsigma " << getDu(index,0,r,sk[r],alpha,1,grid)
//              << " " << getDu(index,1,r,sk[r],alpha,1,grid) << " " 
//              << getDu(index,2,r,sk[r],alpha,1,grid) << endl;
//         cout << "   info = " << alpha << " " << r << " " << sk[r] << " " 
//              << index[0] << " " << index[1] << " " << index[2] << endl;
// form matrix
//         getjumpuxx(jumpu,jumpucoef,jumpuxxcoef,index,r,sk,alpha,thesign,normal,D1, 
//                    D2,S,pb,grid);
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
// form Dn dot product with various vectors
         getDtau(Dtau,index,r,sk[r],alpha,grid);
         getD2tau(D2tau,index,r,sk[r],alpha,grid);
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
// form right hand side vector: b0 without u's 
         for (n = 0; n < grid.dim; n++)
            b0[n] = (sigma/ethere-Dtaudot[0])*dotDndot[n]+dotD2taudot[n];
         b0[grid.dim] = getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                        Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2];
         b0[grid.dim+1] = getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                          Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1];
         b0[grid.dim+2] = -jumpfe;
         
         for (n = 0; n < grid.dim; n++)
            realjumpuxx[n] = getD2u(index,n,n,r,sk[r],alpha,1,grid)-
                             getD2u(index,n,n,r,sk[r],alpha,-1,grid);
         realjumpuxx[grid.dim] = getD2u(index,0,1,r,sk[r],alpha,1,grid)-
                                 getD2u(index,0,1,r,sk[r],alpha,-1,grid);
         realjumpuxx[grid.dim+1] = getD2u(index,1,2,r,sk[r],alpha,1,grid)-
                                   getD2u(index,1,2,r,sk[r],alpha,-1,grid);
         realjumpuxx[grid.dim+2] = getD2u(index,0,2,r,sk[r],alpha,1,grid)-
                                   getD2u(index,0,2,r,sk[r],alpha,-1,grid);
         for (n = 0; n < grid.dim; n++)
            realux[n] = getDu(index,n,r,sk[r],alpha,thesign,grid);
         dotrealuxxdot[0] = 0.0;
         dotrealuxxdot[1] = 0.0;
         for (n = 0; n < grid.dim; n++)
            for (m = 0; m < grid.dim; m++)
            {
               dotrealuxxdot[0] += normal[n]*
                                   getD2u(index,n,m,r,sk[r],alpha,thesign,grid)*
                                   tangent1[m];
               dotrealuxxdot[1] += normal[n]*
                                   getD2u(index,n,m,r,sk[r],alpha,thesign,grid)*
                                   tangent2[m];
            }
         cout << "   checking " << getdotprod(B[0],realjumpuxx,2*grid.dim) << " " 
              << (sigma/ethere+
                  thesign*(ethere-ehere)/ethere*getdotprod(realux,normal,grid.dim)-
                  Dtaudot[0])*dotDndot[0]+dotD2taudot[0] << endl;
         cout << "   checking " << getdotprod(B[1],realjumpuxx,2*grid.dim) << " " 
              << (sigma/ethere+
                  thesign*(ethere-ehere)/ethere*getdotprod(realux,normal,grid.dim)-
                  Dtaudot[0])*dotDndot[1]+dotD2taudot[1] << endl;
         cout << "   checking " << getdotprod(B[2],realjumpuxx,2*grid.dim) << " " 
              << (sigma/ethere+
                  thesign*(ethere-ehere)/ethere*getdotprod(realux,normal,grid.dim)-
                  Dtaudot[0])*dotDndot[2]+dotD2taudot[2] << endl;
         cout << "   checking " << getdotprod(B[3],realjumpuxx,2*grid.dim) << " " 
              << (getdotprod(Dsigma,tangent1,grid.dim)+
                  thesign*(ethere-ehere)*(dotrealuxxdot[0]+
                                          getdotprod(realux,Dndott,grid.dim)))/ethere-
                  Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2] << endl;
         cout << "   checking " << getdotprod(B[4],realjumpuxx,2*grid.dim) << " " 
              << (getdotprod(Dsigma,tangent2,grid.dim)+
                  thesign*(ethere-ehere)*(dotrealuxxdot[1]+
                                          getdotprod(realux,Dndots,grid.dim)))/ethere-
                  Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1] << endl;
         cout << "   checking " << getdotprod(B[5],realjumpuxx,2*grid.dim) << " " 
              << -jumpfe << endl;

// form right hand side vector: bcoef with u's
         for (j = 0; j < grid.dim; j++)
         {
            for (s = 0; s < grid.dim; s++)
               uxxcoef[j][s] = 0.0;
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] < 3)
            {
               setvalarray(ucoef[j],sindex,0.0);

               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 1;
         }
         getDu(ucoef,uxxcoef,index,r,sk[r],alpha,thesign,sk,D1,D2,grid);
         getDu(approxux,index,r,sk[r],alpha,thesign,sk,D1,D2,grid);
         cout << approxux[0] << " " << approxux[1] << " " << approxux[2] << endl;
         cout << evalcoef(0.0,ucoef[0],uxxcoef[0],index,S,grid) << " "
              << evalcoef(0.0,ucoef[1],uxxcoef[1],index,S,grid) << " "
              << evalcoef(0.0,ucoef[2],uxxcoef[2],index,S,grid) << endl;
         exit(1);
         cout << "compare to checking " 
              << (sigma/ethere+
                  thesign*(ethere-ehere)/ethere*getdotprod(approxux,normal,grid.dim)-
                  Dtaudot[0])*dotDndot[1]+dotD2taudot[1] << " "
              << (sigma/ethere-Dtaudot[0])*dotDndot[1]+dotD2taudot[1] << endl;
         realvalue = 0.0;
         for (j = 0; j < grid.dim; j++)
            approxux[j] = 0.0;
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            for (j = 0; j < grid.dim; j++)
               tindex[j] = index[j]-1+sindex[j];
            for (n = 0; n < grid.dim; n++)
               for (m = 0; m < grid.dim; m++)
               {
                  if (theorder == 2)
                     setvalarray(bcoef[n],sindex,evalarray(bcoef[n],sindex)+
                                                 thesign*(ethere-ehere)/ethere*
                                                 dotDndot[n]*normal[m]*
                                                 evalarray(ucoef[m],sindex));
                  else
                     setvalarray(bcoef[n],sindex,evalarray(bcoef[n],sindex)+
                                                 thesign*(ethere-ehere)/ethere*
                                                 dotDndot[n]*normal[m]*
                                                 evalarray(D1[m],sindex));
//                  if (n == 0)
//                     realvalue += thesign*(ethere-ehere)/ethere*dotDndot[n]*
//                                  normal[m]*evalarray(D1[m],sindex)*
//                                  getu(tindex,0,0,0.0,thesign,grid);
                  if (n == 0)
                     approxux[m] += evalarray(D1[m],sindex)*
                                    getu(tindex,0,0,0.0,thesign,grid);
               }
            for (m = 0; m < grid.dim; m++)
            {
               if (theorder == 2)
               {
                  setvalarray(bcoef[grid.dim],sindex,
                              evalarray(bcoef[grid.dim],sindex)+
                              thesign*(ethere-ehere)/ethere*Dndott[m]*
                              evalarray(ucoef[m],sindex));
                  setvalarray(bcoef[grid.dim+1],sindex,
                              evalarray(bcoef[grid.dim+1],sindex)+
                              thesign*(ethere-ehere)/ethere*Dndots[m]*
                              evalarray(ucoef[m],sindex));
               }
               else
               {
                  setvalarray(bcoef[grid.dim],sindex,
                              evalarray(bcoef[grid.dim],sindex)+
                              thesign*(ethere-ehere)/ethere*Dndott[m]*
                              evalarray(D1[m],sindex));
                  setvalarray(bcoef[grid.dim+1],sindex,
                              evalarray(bcoef[grid.dim+1],sindex)+
                              thesign*(ethere-ehere)/ethere*Dndots[m]*
                              evalarray(D1[m],sindex));
               }
               for (n = 0; n < grid.dim; n++)
                  if (n != m)
                  {
                     setvalarray(bcoef[grid.dim],sindex,
                                 evalarray(bcoef[grid.dim],sindex)+
                                 thesign*(ethere-ehere)/ethere*
                                 normal[m]*tangent1[n]*
                                 evalarray(D2[min(m,n)][max(m,n)],sindex));
                     setvalarray(bcoef[grid.dim+1],sindex,
                                 evalarray(bcoef[grid.dim+1],sindex)+
                                 thesign*(ethere-ehere)/ethere*
                                 normal[m]*tangent2[n]*
                                 evalarray(D2[min(m,n)][max(m,n)],sindex));
                     realvalue += thesign*(ethere-ehere)/ethere*
                                  normal[m]*tangent1[n]*
                                  evalarray(D2[min(m,n)][max(m,n)],sindex)*
                                  getu(tindex,0,0,0.0,thesign,grid);
//                     realvalue += thesign*(ethere-ehere)/ethere*
//                                  normal[m]*tangent2[n]*
//                                  evalarray(D2[min(m,n)][max(m,n)],sindex);
                  }
            }
            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;

         cout << "   alpha = " << alpha << endl;
         cout << "   realvalue = " << realvalue << " " 
              << thesign*(ethere-ehere)/ethere*dotrealuxxdot[0] << endl;
//              << getdotprod(realux,normal,grid.dim) << endl;
//                << thesign*(ethere-ehere)/ethere*
//                   getdotprod(realux,normal,grid.dim)*
//                   dotDndot[0] << endl;
         cout << "   realvalue = " << b0[0] << " "
              << (sigma/ethere-Dtaudot[0])*dotDndot[0]+dotD2taudot[0] << endl;
         cout << "   realvalue = " << approxux[0] << " " << approxux[1] << " " 
              << approxux[2] << endl;
         cout << "   realvalue = " << getDu(index,0,r,sk[r],alpha,thesign,grid) << " " 
              << getDu(index,1,r,sk[r],alpha,thesign,grid) << " " 
              << getDu(index,2,r,sk[r],alpha,thesign,grid) << endl;
         for (j = 0; j < grid.dim; j++)
            tindex[j] = index[j];
         tindex[r] -= sk[r];
         cout << "   realvalue = " << (getu(index,0,0,0.0,thesign,grid)-
                                       getu(tindex,0,0,0.0,thesign,grid))/grid.dx[r]
              << endl;
         cout << "   realvalue = " << approxux[0]*normal[0]+approxux[1]*normal[1]+
                                      approxux[2]*normal[2] << " "
              << getdotprod(realux,normal,grid.dim) << endl;
         cout << "   realvalue = " << thesign*(ethere-ehere)/ethere*dotDndot[0]*
                                      (approxux[0]*normal[0]+approxux[1]*normal[1]+
                                       approxux[2]*normal[2]) << " "
              << thesign*(ethere-ehere)/ethere*dotDndot[0]*
                 getdotprod(realux,normal,grid.dim) << endl;
         cout << "   realvalue = " << thesign*(ethere-ehere)/ethere*dotDndot[1]*
                                      (approxux[0]*normal[0]+approxux[1]*normal[1]+
                                       approxux[2]*normal[2]) << " "
              << thesign*(ethere-ehere)/ethere*dotDndot[1]*
                 getdotprod(realux,normal,grid.dim) << endl;
         cout << "   realvalue = " << thesign*(ethere-ehere)/ethere*dotDndot[2]*
                                      (approxux[0]*normal[0]+approxux[1]*normal[1]+
                                       approxux[2]*normal[2]) << " "
              << thesign*(ethere-ehere)/ethere*dotDndot[2]*
                 getdotprod(realux,normal,grid.dim) << endl;

         for (j = 0; j < 2*grid.dim; j++)
            realb[j] = b0[j];
         realvalue = realb[1];
         cout << realb[1] << " " << realvalue << endl;
         for (j = 0; j < grid.dim; j++)
            jumpuxx[j] = 0.0;
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            for (j = 0; j < grid.dim; j++)
               tindex[j] = index[j]-1+sindex[j];
            for (j = 0; j < 2*grid.dim; j++)
               realb[j] += evalarray(bcoef[j],sindex)*getu(tindex,0,0,0.0,thesign,grid);
            
            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         cout << realb[1] << " " << realvalue << endl;
         gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;
         cout << "uxxcoef = " << endl;
         for (m = 0; m < grid.dim; m++)
         {
            for (n = 0; n < grid.dim; n++)
               cout << uxxcoef[m][n] << " ";
            cout << endl;
         }
         for (m = 0; m < grid.dim; m++)
         {
            if (theorder == 2)
            {
               for (n = 0; n < grid.dim; n++)
                  c0[n] = 0.0;
               c0[grid.dim] = thesign*(ethere-ehere)/ethere*normal[m]*tangent1[m];
               c0[grid.dim+1] = thesign*(ethere-ehere)/ethere*normal[m]*tangent2[m];
               for (s = 0; s < grid.dim; s++)
               {
                  for (n = 0; n < grid.dim; n++)
                     c0[n] += thesign*(ethere-ehere)/ethere*dotDndot[n]*
                              normal[s]*uxxcoef[s][m];
                  c0[grid.dim] += thesign*(ethere-ehere)/ethere*Dndott[s]*uxxcoef[s][m];
                  c0[grid.dim+1] += thesign*(ethere-ehere)/ethere*Dndots[s]*uxxcoef[s][m];
               }
            }
            else
            {
               for (n = 0; n < grid.dim; n++)
                  c0[n] = 0.0;
               c0[grid.dim] = thesign*(ethere-ehere)/ethere*normal[m]*tangent1[m];
               c0[grid.dim+1] = thesign*(ethere-ehere)/ethere*normal[m]*tangent2[m];
            }
            c0[grid.dim+2] = 0.0;
            for (j = 0; j < 2*grid.dim; j++)
               realb[j] += c0[j]*getD2u(index,m,m,0,0,0.0,thesign,grid);
            for (n = 0; n < 2*grid.dim; n++)
               temp[n] = c0[n];
            forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
            for (j = 0; j < grid.dim; j++)
               jumpuxx[j] += temp[j]*getD2u(index,m,m,0,0,0.0,thesign,grid);
            cout << "   CHECK NEW PROG " << m << " " << temp[r] << " " 
                 << jumpuxxcoef[m] << endl;
         }
         cout << realb[1] << " " << realvalue << endl;
         cout << "   checking ";
         for (j = 0; j < 2*grid.dim; j++)
            cout << realb[j] << " ";
         cout << endl;
               
         for (j = 0; j < 2*grid.dim; j++)
         {
            for (s = 0; s < 2*grid.dim; s++)
               cout << B[j][s] << " ";
            cout << endl;
         }

// form part of d0 and dcoef's rth entry 
//         gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);
         forwardbacksub0(temp,b0,LU,PLR,PLC,2*grid.dim-1);
         for (j = 0; j < grid.dim; j++)
            jumpuxx[j] += temp[j];
         cout << "   CHECK NEW PROG " << temp[r] << " " << jumpu << endl;
         d0[r] = -0.5*beta*beta*temp[r];
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            for (n = 0; n < 2*grid.dim; n++)
               temp[n] = evalarray(bcoef[n],sindex);
            forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
            cout << "   CHECK NEW PROG " << sindex[0] << " " << sindex[1] << " " 
                 << sindex[2] << " " << temp[r] << " " << evalarray(jumpucoef,sindex) 
                 << endl;
            setvalarray(dcoef[r],sindex,-0.5*beta*beta*temp[r]);
            for (j = 0; j < grid.dim; j++)
               tindex[j] = index[j]-1+sindex[j];
            for (s = 0; s < grid.dim; s++)
               jumpuxx[s] += temp[s]*getu(tindex,0,0,0.0,thesign,grid);

            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;
         for (s = 0; s < grid.dim; s++)
            tindex[s] = index[s];
         cout << "jumpuxx = " << realjumpuxx[0] << " " << realjumpuxx[1] << " " 
              << realjumpuxx[2] << endl;
         cout << "jumpuxx = " << jumpuxx[0] << " " << jumpuxx[1] << " " 
              << jumpuxx[2] << endl;

// form matrix's rth row
         for (m = 0; m < grid.dim; m++)
         {
            if (theorder == 2)
            {
               for (n = 0; n < grid.dim; n++)
                  c0[n] = 0.0;
               c0[grid.dim] = thesign*(ethere-ehere)/ethere*normal[m]*tangent1[m];
               c0[grid.dim+1] = thesign*(ethere-ehere)/ethere*normal[m]*tangent2[m];
               for (s = 0; s < grid.dim; s++)
               {
                  for (n = 0; n < grid.dim; n++)
                     c0[n] += thesign*(ethere-ehere)/ethere*dotDndot[n]*
                              normal[s]*uxxcoef[s][m];
                  c0[grid.dim] += thesign*(ethere-ehere)/ethere*Dndott[s]*uxxcoef[s][m];
                  c0[grid.dim+1] += thesign*(ethere-ehere)/ethere*Dndots[s]*uxxcoef[s][m];
               }
            }
            else
            {
               for (n = 0; n < grid.dim; n++)
                  c0[n] = 0.0;
               c0[grid.dim] = thesign*(ethere-ehere)/ethere*normal[m]*tangent1[m];
               c0[grid.dim+1] = thesign*(ethere-ehere)/ethere*normal[m]*tangent2[m];
            }
            c0[grid.dim+2] = 0.0;
            forwardbacksub0(temp,c0,LU,PLR,PLC,2*grid.dim-1);
            if (m != r)
               G[r][m] = -0.5*sk[r]*sk[m]*beta*(ethere-ehere)/ethere*normal[r]*normal[m]+
                         0.5*beta*beta*temp[r];
            else
               G[r][r] = -alpha*beta*thesign*(ethere-ehere)/ethere*normal[r]*normal[r]+
                         (alpha+0.5)*(alpha+beta)-
                         0.5*beta*(ethere-ehere)/ethere*normal[r]*normal[r]+
                         0.5*beta*beta*temp[r];
         }

// form rest of vector's rth entry
         gettau(tau,index,r,sk[r],alpha,grid);
         cout << "tau = " << tau << endl;
         d0[r] += thesign*sk[r]*beta*
                  (sigma/ethere*normal[r]+
                   Dtaudot[1]*tangent1[r]+Dtaudot[2]*tangent2[r])/grid.dx[r]+
                  thesign*tau/(grid.dx[r]*grid.dx[r]);
         sindex[r] = 1+sk[r];
         for (j = 0; j < grid.dim; j++)
            tindex[j] = index[j]-1+sindex[j];
         setvalarray(dcoef[r],sindex,evalarray(dcoef[r],sindex)+
                                    1.0/(grid.dx[r]*grid.dx[r]));
         sindex[r] = 1;
         setvalarray(dcoef[r],sindex,evalarray(dcoef[r],sindex)-
                                    1.0/(grid.dx[r]*grid.dx[r]));

         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            for (m = 0; m < grid.dim; m++)
            {
               setvalarray(dcoef[r],sindex,evalarray(dcoef[r],sindex)+
                                           sk[r]*beta*(ethere-ehere)/ethere*
                                           normal[r]*normal[m]*evalarray(D1[m],sindex)/
                                           grid.dx[r]);
               if (m != r)
                  setvalarray(dcoef[r],sindex,evalarray(dcoef[r],sindex)+
                                              alpha*beta*(ethere-ehere)/ethere*
                                              normal[r]*normal[m]*
                                              evalarray(D2[min(r,m)][max(r,m)],sindex));
            }
            setvalarray(dcoef[r],sindex,evalarray(dcoef[r],sindex)-
                                        sk[r]*(alpha+beta)*evalarray(D1[r],sindex)/
                                        grid.dx[r]);
            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;

         for (j = 0; j < grid.dim; j++)
            tindex[j] = index[j];
         tindex[r] = index[r]+sk[r];
/*
         reald[r] = -thesign*(tau+
                              sk[r]*beta*grid.dx[r]*
                              (getDu(index,r,r,sk[r],alpha,1,grid)-
                               getDu(index,r,r,sk[r],alpha,-1,grid))+
                              0.5*grid.dx[r]*grid.dx[r]*beta*beta*
                              (getD2u(index,r,r,r,sk[r],alpha,1,grid)-
                               getD2u(index,r,r,r,sk[r],alpha,-1,grid)))+
                     sk[r]*(alpha+beta)*grid.dx[r]*
                     getDu(index,r,r,sk[r],alpha,thesign,grid)+
                     0.5*grid.dx[r]*grid.dx[r]*(beta*beta-alpha*alpha)*
                     getD2u(index,r,r,r,sk[r],alpha,thesign,grid)-
                     (getu(tindex,0,0,0.0,-thesign,grid)-
                      getu(index,0,0,0.0,thesign,grid));
*/
         getDu(approxux,index,r,sk[r],alpha,thesign,sk,D1,D2,grid);
         cout << getDu(index,r,r,sk[r],alpha,1,grid)-
                 getDu(index,r,r,sk[r],alpha,-1,grid) << " "
              << (sigma/ethere+thesign*(ethere-ehere)/ethere*
                               (getDu(index,0,r,sk[r],alpha,thesign,grid)*normal[0]+
                                getDu(index,1,r,sk[r],alpha,thesign,grid)*normal[1]+
                                getDu(index,2,r,sk[r],alpha,thesign,grid)*normal[2]))*
                              normal[r]+Dtaudot[1]*tangent1[r] << " "
              << (sigma/ethere+thesign*(ethere-ehere)/ethere*
                               (approxux[0]*normal[0]+approxux[1]*normal[1]+
                                approxux[2]*normal[2]))*
                              normal[r]+Dtaudot[1]*tangent1[r] << endl;
/*
         reald[r] = -thesign*(tau+
                              sk[r]*beta*grid.dx[r]*
                              ((sigma/ethere+thesign*(ethere-ehere)/ethere*
                                (getDu(index,0,r,sk[r],alpha,thesign,grid)*normal[0]+
                                 getDu(index,1,r,sk[r],alpha,thesign,grid)*normal[1]+
                                 getDu(index,2,r,sk[r],alpha,thesign,grid)*normal[2]))*
                               normal[r]+Dtaudot[1]*tangent1[r])+
                              0.5*grid.dx[r]*grid.dx[r]*beta*beta*
                              (getD2u(index,r,r,r,sk[r],alpha,1,grid)-
                               getD2u(index,r,r,r,sk[r],alpha,-1,grid)))+
                     sk[r]*(alpha+beta)*grid.dx[r]*
                     getDu(index,r,r,sk[r],alpha,thesign,grid)+
                     0.5*grid.dx[r]*grid.dx[r]*(beta*beta-alpha*alpha)*
                     getD2u(index,r,r,r,sk[r],alpha,thesign,grid)-
                     (getu(tindex,0,0,0.0,-thesign,grid)-
                      getu(index,0,0,0.0,thesign,grid));
*/
         cout << getDu(index,r,r,sk[r],alpha,thesign,grid) << " "
              << approxux[r] << endl;
         cout << getD2u(index,r,r,r,sk[r],alpha,thesign,grid) << " "
              << getD2u(index,r,r,0,0,0.0,thesign,grid) << endl;
         reald[r] = -thesign*(tau+
                              sk[r]*beta*grid.dx[r]*
                              ((sigma/ethere+thesign*(ethere-ehere)/ethere*
                                (approxux[0]*normal[0]+approxux[1]*normal[1]+
                                 approxux[2]*normal[2]))*
                               normal[r]+Dtaudot[1]*tangent1[r])+
                              0.5*grid.dx[r]*grid.dx[r]*beta*beta*
                              jumpuxx[r])+
                     sk[r]*(alpha+beta)*grid.dx[r]*approxux[r]+
                     0.5*grid.dx[r]*grid.dx[r]*(beta*beta-alpha*alpha)*
                     getD2u(index,r,r,0,0,0.0,thesign,grid)-
                     (getu(tindex,0,0,0.0,-thesign,grid)-
                      getu(index,0,0,0.0,thesign,grid));
         cout << "   ERROR = " << reald[r] << endl;
         for (j = 0; j < grid.dim; j++)
            tindex[j] = index[j];
      }
      else
      {
         d0[r] = 0.0;
         for (s = 0; s < grid.dim; s++)
            G[r][s] = 0.0;
         G[r][r] = 1.0;
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;
         setvalarray(dcoef[r],sindex,-2.0/(grid.dx[r]*grid.dx[r]));
         for (s = -1; s <= 1; s += 2)
         {
            sindex[r] = 1+s;
            setvalarray(dcoef[r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
         }   
         sindex[r] = 1;
      }
      rindex[r] = index[r];
   }

   double LHS, RHS;
   for (r = 0; r < grid.dim; r++)
   {
      LHS = 0.0;
      for (s = 0; s <= grid.dim; s++)
         LHS += G[r][s]*getD2u(index,s,s,0,0,0.0,thesign,grid);
      RHS = d0[r];
      cout << "d0 " << d0[r] << endl;
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 0;
      while (sindex[0] < 3)
      {
         for (j = 0; j < grid.dim; j++)
            tindex[j] = index[j]-1+sindex[j];
         if (evalarray(S,tindex) < 0.0)
            RHS += evalarray(dcoef[r],sindex)*getu(tindex,0,0,0.0,-1,grid);
         else
            RHS += evalarray(dcoef[r],sindex)*getu(tindex,0,0,0.0,1,grid);
         cout << "d " << sindex[0] << " " << sindex[1] << " " << sindex[2] << " " 
              << evalarray(dcoef[r],sindex) << endl;
         (sindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
         {
            sindex[i] = 0;
            (sindex[i-1])++;
         }
      }
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 1;
      cout << "LHS = " << LHS << endl;
      cout << "RHS = " << RHS << endl;
   }

   double uxxexact, uxxapprox[grid.dim];

   sub2coord(x,index,grid);
   uxxexact = -1.0;
   if (thesign < 0.0)
      for (r = 0; r < grid.dim; r++)
         uxxexact *= sin(x[r]);
   else
      for (r = 0; r < grid.dim; r++)
         uxxexact *= cos(x[r]);

   cout << "d0 = " << d0[0] << " " << d0[1] << " " << d0[2] << endl;
   gecp0(LU,PLR,PLC,G,grid.dim-1,grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,grid.dim-1);
   value = 0.0;
   for (n = 0; n < grid.dim; n++)
      value += ehere*temp[n];
   setvalarray(b,index,evalarray(b,index)+value);
   for (n = 0; n < grid.dim; n++)
   {
      uxxapprox[n] = temp[n];
//      if (n == 1)
//         cout << temp[n] << endl;
   }

   for (s = 0; s < grid.dim; s++)
      tindex[s] = index[s];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      for (n = 0; n < grid.dim; n++)
         temp[n] = evalarray(dcoef[n],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-1+sindex[s];
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value -= ehere*temp[n];
      sparse2(index,tindex,A,value,grid);
      for (n = 0; n < grid.dim; n++)
      {
         if (evalarray(S,tindex) < 0.0)
            uxxapprox[n] += temp[n]*getu(tindex,0,0,0.0,-1,grid);
         else
            uxxapprox[n] += temp[n]*getu(tindex,0,0,0.0,1,grid);
//         if (n == 1)
//            cout << "   " << temp[n] << " " << getu(tindex,0,0,0.0,thesign,grid) 
//                 << " at " << tindex[0] << " " << tindex[1] << " " << tindex[2] 
//                 << " " << uxxapprox[n] << endl;
      }

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;
   for (s = 0; s < grid.dim; s++)
      tindex[s] = index[s];

   SparseElt2 *current;
   current = evalarray(A,index);
//   outputsparserow2(current,grid);
   for (r = 0; r < grid.dim; r++)
   {
      for (s = 0; s < grid.dim; s++)
         cout << G[r][s] << " ";
      cout << endl;
   }
   cout << "uxx = " << uxxexact << " " << uxxapprox[0] << " " << uxxapprox[1] << " "
        << uxxapprox[2] << endl;
   exit(1);

   free_matrix(G,grid.dim-1,grid.dim-1);
   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      free_matrix(D1[r],2,2,2);
      for (s = 0; s < grid.dim; s++)
         free_matrix(D2[r][s],2,2,2);
      free_matrix(dcoef[r],2,2,2);
   }
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(bcoef[r],2,2,2);
}
//cim3 is main idea in cim12.pdf, p28, without higer order CIM1 points
void cim3again(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
               PBData &pb, GridData &grid)
{
   int r, s, t, i, m, n;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];
   double d0[grid.dim];
   double ***dcoef[grid.dim];
//   double ***D1[grid.dim], ***D2[grid.dim][grid.dim];
   double ***D1[grid.dim], ***D2[grid.dim][3];
   double **LU, **G, value;
//   int sk2[grid.dim][grid.dim][4];
   int sk2[grid.dim][3][4];
   int PLR[grid.dim], PLC[grid.dim];
   double alpha, beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[grid.dim];
   int sk[grid.dim];
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;
   double jumpD1u, ***jumpD1ucoef, jumpD1uxxcoef[grid.dim];
   double jumpD2u, ***jumpD2ucoef, jumpD2uxxcoef[grid.dim];
   double ***D1ucoef[grid.dim], **D1uxxcoef;
   

//   getchar();
   double zerouxxcoef[grid.dim];
   double ***zeroucoef;
   zeroucoef = matrix(2,2,2);
   for (s = 0; s < grid.dim; s++) 
      zerouxxcoef[s] = 0.0;
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      setvalarray(zeroucoef,sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;
   value = 0.0;

   LU = matrix(grid.dim-1,grid.dim-1);
   G = matrix(grid.dim-1,grid.dim-1);
   jumpD1ucoef = matrix(2,2,2);
   jumpD2ucoef = matrix(2,2,2);
   for (r = 0; r < grid.dim; r++)
   {
      D1[r] = matrix(2,2,2);
      for (s = 0; s < grid.dim; s++)
         D2[r][s] = matrix(2,2,2);
      dcoef[r] = matrix(2,2,2);
   }
   for (r = 0; r < grid.dim; r++)
      D1ucoef[r] = matrix(2,2,2);
   D1uxxcoef = matrix(grid.dim-1,grid.dim-1);

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
      sk[r] = gamma[r][1]-gamma[r][0];
// get sk2
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
         getsk2(sk2[m][n],m,n,index,S,grid);

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;
   for (r = 0; r < grid.dim; r++)
      if (abs(sk[r]) == 1)
      {
         rindex[r] = index[r]+sk[r];
// getting first and second derivatives,
         getD1(D1,sk,grid);
         getD2(D2,sk2,grid);
// getting interface info
         getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
         beta = 1.0-alpha;
// get Du 2nd order
         getDu(D1ucoef,D1uxxcoef,index,r,sk[r],alpha,thesign,sk,D1,D2,grid);
// getting jump of uxx in r
         gettau(tau,index,r,sk[r],alpha,grid);
         getjumpux(jumpD1u,jumpD1ucoef,jumpD1uxxcoef,index,r,sk,alpha,thesign,normal,
                   tangent,D1ucoef,D1uxxcoef,S,pb,grid);
         if (globorder == 1)
            getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,index,r,sk,alpha,thesign,
                       normal,D1,D2,S,pb,grid);
         else
            getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,index,r,sk,alpha,thesign,
                       normal,D1ucoef,D1uxxcoef,D2,S,pb,grid);
         if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
         {
            cout << "rstar = " << r << " sstar = " << sk[r] << endl;
            cout << "alpha = " << alpha << endl;
            cout << "sk = " << sk[0] << " " << sk[1] << " " << sk[2] << endl;
            cout << "   check jumpux " 
                 << evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxxcoef,index,S,grid) << " "
                 << getDu(index,r,r,sk[r],alpha,1,grid)-
                    getDu(index,r,r,sk[r],alpha,-1,grid) << endl;
            cout << "   check jumpuxx " 
                 << evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,index,S,grid) << " "
                 << getD2u(index,r,r,r,sk[r],alpha,1,grid)-
                    getD2u(index,r,r,r,sk[r],alpha,-1,grid) << endl;
            cout << "   check ux " 
                 << evalcoef(0.0,D1ucoef[0],D1uxxcoef[0],index,S,grid) << " "
                 << getDu(index,0,r,sk[r],alpha,thesign,grid) << " "
                 << evalcoef(0.0,D1[0],zerouxxcoef,index,S,grid) << endl;
            cout << "   check ux " 
                 << evalcoef(0.0,D1ucoef[1],D1uxxcoef[1],index,S,grid) << " "
                 << getDu(index,1,r,sk[r],alpha,thesign,grid) << " "
                 << evalcoef(0.0,D1[1],zerouxxcoef,index,S,grid) << endl;
            cout << "   check ux " 
                 << evalcoef(0.0,D1ucoef[2],D1uxxcoef[2],index,S,grid) << " "
                 << getDu(index,2,r,sk[r],alpha,thesign,grid) << " "
                 << evalcoef(0.0,D1[2],zerouxxcoef,index,S,grid) << endl;
            value = getu(index,r,sk[r],1.0,-thesign,grid)-
                    getu(index,0,0,0.0,thesign,grid)+
                    thesign*(tau+
                             sk[r]*beta*grid.dx[r]*
                             evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxxcoef,index,S,grid)+
                             0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                             evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,index,S,grid))-
                    sk[r]*grid.dx[r]*(alpha+beta)*
                    evalcoef(0.0,D1ucoef[r],D1uxxcoef[r],index,S,grid)-
                    0.5*grid.dx[r]*grid.dx[r]*(beta*beta-alpha*alpha)*
                    getD2u(index,r,r,0,0,0.0,thesign,grid);
            cout << value << " " << value/(grid.dx[r]*grid.dx[r]) << endl;
            value = thesign*(sk[r]*beta*grid.dx[r]*
                             evalcoef(0.0,zeroucoef,jumpD1uxxcoef,index,S,grid)+
                             0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                             evalcoef(0.0,zeroucoef,jumpD2uxxcoef,index,S,grid))-
                    sk[r]*grid.dx[r]*(alpha+beta)*
                    evalcoef(0.0,zeroucoef,D1uxxcoef[r],index,S,grid)-
                    0.5*grid.dx[r]*grid.dx[r]*(beta*beta-alpha*alpha)*
                    getD2u(index,r,r,0,0,0.0,thesign,grid);
            cout << value << " " << value/(grid.dx[r]*grid.dx[r]) << endl;
         }
// form d0 and dcoef's rth entry. G [uxx uyy uzz] = d0 + dcoeff[r][3x3x3] uij, cim12.pdf p25
         d0[r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+sk[r]*beta*jumpD1u/grid.dx[r]+
                          0.5*beta*beta*jumpD2u);

         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            setvalarray(dcoef[r],sindex,
                        thesign*(sk[r]*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                 0.5*beta*beta*evalarray(jumpD2ucoef,sindex))-
                        sk[r]*(alpha+beta)*evalarray(D1ucoef[r],sindex)/grid.dx[r]);

            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;
         setvalarray(dcoef[r],sindex,
                     evalarray(dcoef[r],sindex)-1.0/(grid.dx[r]*grid.dx[r]));
         sindex[r] = 1+sk[r];
         setvalarray(dcoef[r],sindex,
                     evalarray(dcoef[r],sindex)+1.0/(grid.dx[r]*grid.dx[r]));
         sindex[r] = 1;
// form matrix's rth row
         for (m = 0; m < grid.dim; m++)
            G[r][m] = -thesign*(sk[r]*beta*jumpD1uxxcoef[m]/grid.dx[r]+
                                0.5*beta*beta*jumpD2uxxcoef[m])+
                      sk[r]*(alpha+beta)*D1uxxcoef[r][m]/grid.dx[r];
         G[r][r] += 0.5*(beta*beta-alpha*alpha);

         rindex[r] = index[r];
         if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
         {
            cout << G[r][0]*getD2u(index,0,0,0,0,0.0,thesign,grid)+
                    G[r][1]*getD2u(index,1,1,0,0,0.0,thesign,grid)+
                    G[r][2]*getD2u(index,2,2,0,0,0.0,thesign,grid) << endl;
            cout << "   check RHS " 
                 << evalcoef(d0[r],dcoef[r],zerouxxcoef,index,S,grid) << endl;
         }
      }
      else
      {
// if interior point
         d0[r] = 0.0;
         for (s = 0; s < grid.dim; s++)
            G[r][s] = 0.0;
         G[r][r] = 1.0;
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            setvalarray(dcoef[r],sindex,0.0);

            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;
         setvalarray(dcoef[r],sindex,-2.0/(grid.dx[r]*grid.dx[r]));
         for (s = -1; s <= 1; s += 2)
         {
            sindex[r] = 1+s;
            setvalarray(dcoef[r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
         }
         sindex[r] = 1;
      }
   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
   {
      for (m = 0; m < grid.dim; m++)
      {
         for (s = 0; s < grid.dim; s++)
            cout << G[m][s] << " ";
         cout << endl;
      }
   }
//   for (m = 0; m < grid.dim; m++)
//   {
//      for (n = 0; n < grid.dim; n++)
//         printf("%4.16f ",G[m][n]);
//      printf("%4.16f\n",evalcoef(d0[m],dcoef[m],zerouxxcoef,index,S,grid));
//   }

// solve for uxx and put in matrix
   int j;
   double uxx[grid.dim];
   for (j = 0; j < grid.dim; j++)
      uxx[j] = 0.0;

   gecp0(LU,PLR,PLC,G,grid.dim-1,grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,grid.dim-1);
   value = 0.0;
   for (n = 0; n < grid.dim; n++)
      value += ehere*temp[n];
//   setvalarray(b,index,evalarray(b,index)+value);
   setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid)+value);
   for (j = 0; j < grid.dim; j++)
      uxx[j] = temp[j];//uxx is calculated by plug in exact u for purpose of checking 

   for (s = 0; s < grid.dim; s++)
      tindex[s] = index[s];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3) //sindex is local index of nbr points
   {
      for (n = 0; n < grid.dim; n++)
         temp[n] = evalarray(dcoef[n],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-1+sindex[s];//tindex is global index of nbr points, 
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value -= ehere*temp[n];
      if (value != 0.0)
         sparse2(index,tindex,A,value,grid);
      for (j = 0; j < grid.dim; j++)
         if (evalarray(S,tindex) < 0.0)
            uxx[j] += temp[j]*getu(tindex,0,0,0.0,-1,grid);
         else
            uxx[j] += temp[j]*getu(tindex,0,0,0.0,1,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
   {
      cout << "index = " << index[0] << " " << index[1] << " " << index[2] << endl;
      cout << "uxx = " << uxx[0] << " " << uxx[1] << " " << uxx[2] << endl;
      cout << "uxx = " << getD2u(index,0,0,0,0,0.0,thesign,grid) << " " 
           << getD2u(index,1,1,0,0,0.0,thesign,grid) << " "
           << getD2u(index,2,2,0,0,0.0,thesign,grid) << endl;
      cout << "equation = " << ehere*(uxx[0]+uxx[1]+uxx[2])+
                               getf(index,0,0,0.0,thesign,pb,grid) << endl;
   }
   if (globerr < fabs(ehere*(uxx[0]+uxx[1]+uxx[2])+
                      getf(index,0,0,0.0,thesign,pb,grid)))
   {
      globerr = fabs(ehere*(uxx[0]+uxx[1]+uxx[2])+
                     getf(index,0,0,0.0,thesign,pb,grid));
      for (s = 0; s < grid.dim; s++)
         globerrvec3[s] = index[s];
   }

   free_matrix(zeroucoef,2,2,2);
   free_matrix(LU,grid.dim-1,grid.dim-1);
   free_matrix(G,grid.dim-1,grid.dim-1);
   free_matrix(jumpD1ucoef,2,2,2);
   free_matrix(jumpD2ucoef,2,2,2);
   for (r = 0; r < grid.dim; r++)
   {
      free_matrix(D1[r],2,2,2);
      for (s = 0; s < grid.dim; s++)
         free_matrix(D2[r][s],2,2,2);
      free_matrix(dcoef[r],2,2,2);
   }
   for (r = 0; r < grid.dim; r++)
      free_matrix(D1ucoef[r],2,2,2);
   free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1);
   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
      exit(1);
//   getchar();
}
//cim4 improves cim3 by bringing cross-derivative from the other side, cim12.pdf p 36
void cim4(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
          PBData &pb, GridData &grid)
{
   int r, s, t, i, m, n, sk;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];
   double d0[2*grid.dim];
   double ***dcoef[2*grid.dim];
//   double ***D2[grid.dim][grid.dim], D2jumpuxxcoef[grid.dim][grid.dim];
   double ***D2[grid.dim][3], D2jumpuxxcoef[grid.dim][3];
   double **LU, **G, value;
//   int sk2[grid.dim][grid.dim][4];
   int sk2[grid.dim][3][4];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double alpha, beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[2*grid.dim], **temp2;
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;
   double jumpD1u = 0.0, ***jumpD1ucoef, jumpD1uxcoef[grid.dim], **jumpD1uxxcoef;
   double **jumpD2u, *****jumpD2ucoef, ***jumpD2uxcoef, ***jumpD2uxxcoef;
   double D1u[grid.dim], ***D1ucoef[grid.dim], **D1uxcoef, ***D1uxxcoef;

//   getchar();
   double zerouxxcoef[grid.dim], zerouxcoef[grid.dim];
   double ***zeroucoef;
   zeroucoef = matrix(4,4,4);
   D1uxcoef = matrix(grid.dim-1,grid.dim-1);
   D1uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   jumpD1uxxcoef = matrix(grid.dim-1,grid.dim-1);
   temp2 = matrix(grid.dim-1,grid.dim-1);
   for (s = 0; s < grid.dim; s++) 
   {
      zerouxcoef[s] = 0.0;
      zerouxxcoef[s] = 0.0;
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 5)
   {
      setvalarray(zeroucoef,sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 2;
   value = 0.0;

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   G = matrix(2*grid.dim-1,2*grid.dim-1);
   jumpD1ucoef = matrix(4,4,4);
   jumpD2u = matrix(grid.dim-1,grid.dim-1);
   jumpD2uxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   jumpD2uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   jumpD2ucoef = new double ****[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      D1ucoef[r] = matrix(4,4,4);
      jumpD2ucoef[r] = new double ***[grid.dim];
      for (s = 0; s < grid.dim; s++)
      {
         jumpD2ucoef[r][s] = matrix(4,4,4);
         D2[r][s] = matrix(4,4,4);
         D2jumpuxxcoef[r][s] = 0.0;
         temp2[r][s] = 0.0;
      }
   }
   for (r = 0; r < 2*grid.dim; r++)
      dcoef[r] = matrix(4,4,4);

   double ux[grid.dim];
   for (r = 0; r < grid.dim; r++)
      ux[r] = 0.0;

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
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 2;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
            for (s = 0; s < grid.dim; s++)
               D1u[s] = 0.0;
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] < 5)
            {
               for (s = 0; s < grid.dim; s++)
                  setvalarray(D1ucoef[s],sindex,0.0);
               setvalarray(jumpD1ucoef,sindex,0.0);

               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 2;
            
            rindex[r] = index[r]+sk;
// getting interface info
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
               cout << "S = " << evalarray(S,index) << endl;
            getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
            beta = 1.0-alpha;
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
               printf("alpha = %4.16f\n",alpha);
// get sk2 and second derivatives
            for (m = 0; m < grid.dim; m++)
               for (n = m+1; n < grid.dim; n++)
               {//if can get cross-deri at index. should not get in here, already prevented by status classification
                  if (getsk2(sk2[m][n],m,n,index,S,grid)) 
                  {
                     getD2(D2[m][n],m,n,sk2,sindex,grid);
                     if (index[0] == eindex[0] && index[1] == eindex[1] && 
                         index[2] == eindex[2])
                     {
                        cout << "uxx regular" << endl; 
                        cout << evalcoef(0.0,D2[m][n],zerouxcoef,zerouxxcoef,temp2,
                                         index,r,sk,alpha,S,grid) << " " 
                             << getD2u(index,m,n,0,0,0.0,thesign,grid) << " "
                             << getD2u(index,m,n,r,sk,alpha,thesign,grid) << endl;
                     }
                  }
                  else
                  {
//                     getsk2(sk2[m][n],m,n,rindex,S,grid);
                     if (!getsk2(sk2[m][n],m,n,rindex,S,grid)) //if cannot get cross-deri at index and across interface
                     {
                        cout << "found a non cim4 point" << endl;
                        exit(1);
                     }
                     getD2(D2[m][n],D2jumpuxxcoef[m][n],m,n,r,sk,thesign,sk2,grid);
                     if (index[0] == eindex[0] && index[1] == eindex[1] && 
                         index[2] == eindex[2])
                     {
                        cout << "uxx irregular" << endl; 
                        temp2[m][n] = D2jumpuxxcoef[m][n];
                        cout << evalcoef(0.0,D2[m][n],zerouxcoef,zerouxxcoef,temp2,
                                         index,r,sk,alpha,S,grid) << " " 
                             << getD2u(index,m,n,0,0,0.0,thesign,grid) << " "
                             << getD2u(index,m,n,r,sk,alpha,thesign,grid) << endl;
                        temp2[m][n] = 0.0;
                        cout << evalcoef(0.0,D2[m][n],zerouxcoef,zerouxxcoef,temp2,
                                         index,r,sk,alpha,S,grid) << " " 
                             << getD2u(index,m,n,r,sk,1.0,-thesign,grid) << " "
                             << getD2u(rindex,m,n,0,0,0.0,-thesign,grid) << " "
                             << getD2u(index,m,n,r,sk,alpha,-thesign,grid) << endl;
                        cout << getstatus4(S,index,grid) << " " 
                             << getstatus4(S,rindex,grid) << endl;
                     }
                  }
               }
            getDu(D1uxcoef,D1uxxcoef,index,r,sk,alpha,thesign,grid);
            value = 0.0;
            for (s = 0; s < grid.dim; s++)
               value += D1uxcoef[r][s]*getDu(index,s,0,0,0.0,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  value += D1uxxcoef[r][m][n]*getD2u(index,m,n,0,0,0.0,thesign,grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               cout << "ux" << endl;
               cout << value << " " << getDu(index,r,r,sk,alpha,thesign,grid) << endl;
               ux[r] = sk*(getu(index,0,0,0.0,thesign,grid)-
                           getu(index,r,-sk,1.0,thesign,grid))/grid.dx[r]+
                       sk*grid.dx[r]/2.0*getD2u(index,r,r,0,0,0.0,thesign,grid)+
                       sk*grid.dx[r]*alpha*getD2u(index,r,r,0,0,0.0,thesign,grid);
               for (m = 0; m < grid.dim; m++)
                  if (m != r)
                     ux[m] = (getu(index,m,1,1.0,thesign,grid)-
                              getu(index,m,-1,1.0,thesign,grid))/(2.0*grid.dx[r])+
                             sk*grid.dx[r]*alpha*evalcoef(0.0,D2[min(m,r)][max(m,r)],
                                                          zerouxcoef,zerouxxcoef,temp2,
                                                          index,r,sk,alpha,S,grid);
               cout << ux[0] << " " << ux[1] << " " << ux[2] << endl;
               ux[r] = sk*(getu(index,0,0,0.0,thesign,grid)-
                           getu(index,r,-sk,1.0,thesign,grid))/grid.dx[r]+
                       sk*grid.dx[r]/2.0*getD2u(index,r,r,0,0,0.0,thesign,grid);
               for (m = 0; m < grid.dim; m++)
                  if (m != r)
                     ux[m] = (getu(index,m,1,1.0,thesign,grid)-
                              getu(index,m,-1,1.0,thesign,grid))/(2.0*grid.dx[r]);
            }
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha,grid);
            getjumpux(jumpD1u,jumpD1uxcoef,jumpD1uxxcoef,index,r,sk,alpha,thesign,normal,
                      tangent,D1uxcoef,D1uxxcoef,S,pb,grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               value = jumpD1u;
               for (s = 0; s < grid.dim; s++)
                  value += jumpD1uxcoef[s]*getDu(index,s,0,0,0.0,thesign,grid);
               for (m = 0; m < grid.dim; m++)
                  for (n = m; n < grid.dim; n++)
                     value += jumpD1uxxcoef[m][n]*getD2u(index,m,n,0,0,0.0,thesign,grid);
               cout << "jump ux" << endl;
               cout << value << " " << getDu(index,r,r,sk,alpha,1.0,grid)-
                                       getDu(index,r,r,sk,alpha,-1.0,grid) << endl;
               value = jumpD1u;
               for (s = 0; s < grid.dim; s++)
                  value += jumpD1uxcoef[s]*ux[s];
               for (m = 0; m < grid.dim; m++)
                  for (n = m; n < grid.dim; n++)
                     if (n == m)
                        value += jumpD1uxxcoef[m][n]*getD2u(index,m,n,0,0,0.0,thesign,
                                                            grid);
                     else
                        value += jumpD1uxxcoef[m][n]*
                                 evalcoef(0.0,D2[min(m,n)][max(m,n)],zerouxcoef,
                                          zerouxxcoef,temp2,index,r,sk,alpha,S,grid);
               cout << value << endl;
            }
            getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,alpha,
                       thesign,normal,D1uxcoef,D1uxxcoef,D2,D2jumpuxxcoef,S,pb,
                       grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               for (m = 0; m < grid.dim; m++)
                  for (n = m; n < grid.dim; n++)
                  {
                     value = jumpD2u[m][n];
                     for (s = 0; s < grid.dim; s++)
                     {
                        value += jumpD2uxcoef[m][n][s]*getDu(index,s,0,0,0.0,thesign,
                                                             grid);
                        value += jumpD2uxxcoef[m][n][s]*
                                 getD2u(index,s,s,0,0,0.0,thesign,grid);
                     }
                     for (s = 0; s < grid.dim; s++)
                        sindex[s] = 0;
                     while (sindex[0] < 5)
                     {
                        for (s = 0; s < grid.dim; s++)
                           tindex[s] = index[s]-2+sindex[s];
                        if (evalarray(S,tindex) < 0.0)
                           value += evalarray(jumpD2ucoef[m][n],sindex)*
                                    getu(tindex,0,0,0.0,-1,grid);
                        else
                           value += evalarray(jumpD2ucoef[m][n],sindex)*
                                    getu(tindex,0,0,0.0,1,grid);
                        (sindex[grid.dim-1])++;
                        for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
                        {
                           sindex[i] = 0;
                           (sindex[i-1])++;
                        }
                     }
                     for (s = 0; s < grid.dim; s++)
                        sindex[s] = 2;
                     cout << "jump uxx" << endl;
                     cout << value << " " << getD2u(index,m,n,r,sk,alpha,1.0,grid)-
                                             getD2u(index,m,n,r,sk,alpha,-1.0,grid) 
                          << endl;
                     value = jumpD2u[m][n];
                     for (s = 0; s < grid.dim; s++)
                     {
                        value += jumpD2uxcoef[m][n][s]*ux[s];
                        value += jumpD2uxxcoef[m][n][s]*
                                 getD2u(index,s,s,0,0,0.0,thesign,grid);
                     }
                     for (s = 0; s < grid.dim; s++)
                        sindex[s] = 0;
                     while (sindex[0] < 5)
                     {
                        for (s = 0; s < grid.dim; s++)
                           tindex[s] = index[s]-2+sindex[s];
                        if (evalarray(S,tindex) < 0.0)
                           value += evalarray(jumpD2ucoef[m][n],sindex)*
                                    getu(tindex,0,0,0.0,-1,grid);
                        else
                           value += evalarray(jumpD2ucoef[m][n],sindex)*
                                    getu(tindex,0,0,0.0,1,grid);
                        (sindex[grid.dim-1])++;
                        for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
                        {
                           sindex[i] = 0;
                           (sindex[i-1])++;
                        }
                     }
                     for (s = 0; s < grid.dim; s++)
                        sindex[s] = 2;
                     cout << value << endl;
                  }
            }
            recast(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,jumpD2u,jumpD2ucoef,
                   jumpD2uxcoef,jumpD2uxxcoef,D2,D2jumpuxxcoef,grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               value = jumpD1u;
               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 0;
               while (sindex[0] < 5)
               {
                  for (s = 0; s < grid.dim; s++)
                     tindex[s] = index[s]-2+sindex[s];
                  if (evalarray(S,tindex) < 0.0)
                     value += evalarray(jumpD1ucoef,sindex)*
                              getu(tindex,0,0,0.0,-1,grid);
                  else
                     value += evalarray(jumpD1ucoef,sindex)*
                              getu(tindex,0,0,0.0,1,grid);
                  (sindex[grid.dim-1])++;
                  for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
                  {
                     sindex[i] = 0;
                     (sindex[i-1])++;
                  }
               }
               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 2;
               for (s = 0; s < grid.dim; s++)
//               value += jumpD1uxcoef[s]*getDu(index,s,0,0,0.0,thesign,grid);
                  value += jumpD1uxcoef[s]*ux[s];
               for (m = 0; m < grid.dim; m++)
                  for (n = m; n < grid.dim; n++)
                     value += jumpD1uxxcoef[m][n]*getD2u(index,m,n,0,0,0.0,thesign,grid);
               cout << "jump ux" << endl;
               cout << value << " " << getDu(index,r,r,sk,alpha,1.0,grid)-
                                       getDu(index,r,r,sk,alpha,-1.0,grid) << endl;
               value = 0.0;
               for (s = 0; s < grid.dim; s++)
                  value += D1uxcoef[r][s]*getDu(index,s,0,0,0.0,thesign,grid);
               cout << "before recast ux" << endl;
               cout << value << endl;
               for (m = 0; m < grid.dim; m++)
                  for (n = m; n < grid.dim; n++)
                     value += D1uxxcoef[r][m][n]*getD2u(index,m,n,0,0,0.0,thesign,grid);
               cout << value << endl;
            }
            recast(D1u[r],D1ucoef[r],D1uxcoef[r],D1uxxcoef[r],jumpD2u,jumpD2ucoef,
                   jumpD2uxcoef,jumpD2uxxcoef,D2,D2jumpuxxcoef,grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               value = D1u[r];
               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 0;
               while (sindex[0] < 5)
               {
                  for (s = 0; s < grid.dim; s++)
                     tindex[s] = index[s]-2+sindex[s];
                  if (evalarray(S,tindex) < 0.0)
                     value += evalarray(D1ucoef[r],sindex)*
                              getu(tindex,0,0,0.0,-1,grid);
                  else
                     value += evalarray(D1ucoef[r],sindex)*
                              getu(tindex,0,0,0.0,1,grid);
                  (sindex[grid.dim-1])++;
                  for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
                  {
                     sindex[i] = 0;
                     (sindex[i-1])++;
                  }
               }
               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 2;
               for (s = 0; s < grid.dim; s++)
                  value += D1uxcoef[r][s]*getDu(index,s,0,0,0.0,thesign,grid);
               for (m = 0; m < grid.dim; m++)
                  for (n = m; n < grid.dim; n++)
                     value += D1uxxcoef[r][m][n]*getD2u(index,m,n,0,0,0.0,thesign,grid);
               cout << "ux" << endl;
               cout << value << " " << getDu(index,r,r,sk,alpha,thesign,grid) << endl;
            }
// form d0 and dcoef's rth entry. sk=-1, r=1,2,3 for row 123, sk=1, r=1,2,3 for row 456, vector x = [uxx uyy uzz ux uy uz]
            if (globoldwrongcim4)
               d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                                  sk*beta*jumpD1u/grid.dx[r]+
                                                  0.5*beta*beta*jumpD2u[r][r]);
            else
               d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                                  sk*beta*jumpD1u/grid.dx[r]+
                                                  0.5*beta*beta*jumpD2u[r][r])-
                                         sk*D1u[r]/grid.dx[r];

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] < 5)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                           thesign*(sk*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                    0.5*beta*beta*evalarray(jumpD2ucoef[r][r],sindex))-
                           sk*evalarray(D1ucoef[r],sindex)/grid.dx[r]);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 2;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)-
                        1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 2+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)+
                        1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 2;
// form matrix's rth row
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m] = -thesign*
                                           (sk*beta*jumpD1uxxcoef[m][m]/grid.dx[r]+
                                            0.5*beta*beta*jumpD2uxxcoef[r][r][m])+
                                           sk*D1uxxcoef[r][m][m]/grid.dx[r];
            G[grid.dim*(sk+1)/2+r][r] += 0.5*(beta*beta-alpha*alpha);
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m+grid.dim] = -thesign*
                                                    (sk*beta*jumpD1uxcoef[m]/grid.dx[r]+
                                                     0.5*beta*beta*jumpD2uxcoef[r][r][m])+
                                                    sk*D1uxcoef[r][m]/grid.dx[r];
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               value = 0.0;
               for (m = 0; m < grid.dim; m++)
               {
                  value += G[grid.dim*(sk+1)/2+r][m]*
                           getD2u(index,m,m,0,0,0.0,thesign,grid);
                  value += G[grid.dim*(sk+1)/2+r][m+grid.dim]*
                           getDu(index,m,0,0,0.0,thesign,grid);
               }
               for (m = 0; m < grid.dim; m++)
                  for (n = 0; n < grid.dim; n++)
                     temp2[m][n] = 0.0;
               cout << "not interior" << endl;
               cout << "matrix " << value << endl;
               cout << "rhs " << evalcoef(d0[grid.dim*(sk+1)/2+r],
                                          dcoef[grid.dim*(sk+1)/2+r],
                                          zerouxcoef,zerouxxcoef,temp2,index,0,0,0.0,
                                          S,grid) << endl;
            }

            rindex[r] = index[r];
         }
         else
         {
// if interior point, u(r+s) = u(r) + s h u'(r) + 1/2 h^2 u''(r) 
            d0[grid.dim*(sk+1)/2+r] = 0.0;
            for (s = 0; s < 2*grid.dim; s++)
               G[grid.dim*(sk+1)/2+r][s] = 0.0;
            G[grid.dim*(sk+1)/2+r][r] = 0.5;
            G[grid.dim*(sk+1)/2+r][grid.dim+r] = sk/grid.dx[r];
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] < 5)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,0.0);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 2;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,-1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 2+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 2;

            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               value = 0.0;
               for (m = 0; m < grid.dim; m++)
               {
                  value += G[grid.dim*(sk+1)/2+r][m]*
                           getD2u(index,m,m,0,0,0.0,thesign,grid);
                  value += G[grid.dim*(sk+1)/2+r][m+grid.dim]*
                           getDu(index,m,0,0,0.0,thesign,grid);
               }
               cout << "interior" << endl;
               cout << "matrix " << value << endl;
               cout << "rhs " << evalcoef(d0[grid.dim*(sk+1)/2+r],
                                          dcoef[grid.dim*(sk+1)/2+r],
                                          zerouxcoef,zerouxxcoef,temp2,index,0,0,0.0,
                                          S,grid) << endl;
            }
         }

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
   {
      for (m = 0; m < 2*grid.dim; m++)
      {
         for (n = 0; n < 2*grid.dim; n++)
            printf("%4.16f ",G[m][n]);
         printf("%4.16f\n",evalcoef(d0[m],dcoef[m],
                                    zerouxcoef,zerouxxcoef,temp2,index,0,0,0.0,S,grid));
      }
   }
// solve for uxx and put in matrix
   int j;
   double uxx[grid.dim];
   for (j = 0; j < grid.dim; j++)
      uxx[j] = 0.0;

   gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,2*grid.dim-1);
   value = 0.0;
   for (n = 0; n < grid.dim; n++)
      value += ehere*temp[n];
//   setvalarray(b,index,evalarray(b,index)+value);
   setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid)+value);//RHS
   for (j = 0; j < grid.dim; j++)
      uxx[j] = temp[j];
   for (j = 0; j < grid.dim; j++)
      ux[j] = temp[j+grid.dim];

   for (s = 0; s < grid.dim; s++)
      tindex[s] = index[s];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 5)
   {
      for (n = 0; n < 2*grid.dim; n++)
         temp[n] = evalarray(dcoef[n],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);//get [DDu Du] = inv(G) dcoef[1:6][sindex]
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-2+sindex[s];
      value = 0.0;
      for (n = 0; n < grid.dim; n++)// first 3 row is DDu term
         value -= ehere*temp[n];
      if (value != 0.0)
         sparse2(index,tindex,A,value,grid);
       // below for checking
      for (j = 0; j < grid.dim; j++)
         if (evalarray(S,tindex) < 0.0)
            uxx[j] += temp[j]*getu(tindex,0,0,0.0,-1,grid);
         else
            uxx[j] += temp[j]*getu(tindex,0,0,0.0,1,grid);
      for (j = 0; j < grid.dim; j++)
         if (evalarray(S,tindex) < 0.0)
            ux[j] += temp[j+grid.dim]*getu(tindex,0,0,0.0,-1,grid);
         else
            ux[j] += temp[j+grid.dim]*getu(tindex,0,0,0.0,1,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 2;

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
   {
      cout << "cim4" << endl;
      cout << "index = " << index[0] << " " << index[1] << " " << index[2] << endl;
      printf("ux = %4.16f %4.16f %4.16f\n",ux[0],ux[1],ux[2]);
      cout << "ux = " << ux[0] << " " << ux[1] << " " << ux[2] << endl;
      printf("real ux = %4.16f %4.16f %4.16f\n",getDu(index,0,0,0,0.0,thesign,grid),
                                                getDu(index,1,0,0,0.0,thesign,grid),
                                                getDu(index,2,0,0,0.0,thesign,grid));
      cout << "uxx = " << uxx[0] << " " << uxx[1] << " " << uxx[2] << endl;
      printf("real uxx = %4.16f %4.16f %4.16f\n",getD2u(index,0,0,0,0,0.0,thesign,grid),
                                                 getD2u(index,1,1,0,0,0.0,thesign,grid),
                                                 getD2u(index,2,2,0,0,0.0,thesign,grid));
      cout << "real uxx = " << getD2u(index,0,0,0,0,0.0,thesign,grid) << " " 
           << getD2u(index,1,1,0,0,0.0,thesign,grid) << " "
           << getD2u(index,2,2,0,0,0.0,thesign,grid) << endl;
      cout << "equation = " << ehere*(uxx[0]+uxx[1]+uxx[2])+
                               getf(index,0,0,0.0,thesign,pb,grid) << endl;
   }
/*
   if (index[0] == 34 && index[1] == 31 && index[2] == 56)
   {
      cout << ux[0] << " " << ux[1] << " " << ux[2] << endl;
      cout << uxx[0] << " " << uxx[1] << " " << uxx[2] << endl;
      cout << endl; 
      cout << getDu(index,0,0,0,0.0,thesign,grid) << " "
           << getDu(index,1,0,0,0.0,thesign,grid) << " "
           << getDu(index,2,0,0,0.0,thesign,grid) << endl;
      cout << getD2u(index,0,0,0,0,0.0,thesign,grid) << " "
           << getD2u(index,1,1,0,0,0.0,thesign,grid) << " "
           << getD2u(index,2,2,0,0,0.0,thesign,grid) << endl;
      getchar();
   }
*/
   if (globerr < fabs(ehere*(uxx[0]+uxx[1]+uxx[2])+
                      getf(index,0,0,0.0,thesign,pb,grid)))
   {
      globerr = fabs(ehere*(uxx[0]+uxx[1]+uxx[2])+
                     getf(index,0,0,0.0,thesign,pb,grid));
      for (s = 0; s < grid.dim; s++)
         globerrvec3[s] = index[s];
   }

   free_matrix(zeroucoef,4,4,4);
   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(G,2*grid.dim-1,2*grid.dim-1);
   free_matrix(jumpD1ucoef,4,4,4);
   free_matrix(jumpD1uxxcoef,grid.dim-1,grid.dim-1);
   free_matrix(jumpD2u,grid.dim-1,grid.dim-1);
   free_matrix(jumpD2uxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(jumpD2uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D1uxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(temp2,grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      free_matrix(D1ucoef[r],4,4,4);
      for (s = 0; s < grid.dim; s++)
      {
         free_matrix(D2[r][s],4,4,4);
         free_matrix(jumpD2ucoef[r][s],4,4,4);
      }
      delete [] jumpD2ucoef[r];
   }
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(dcoef[r],4,4,4);
   delete [] jumpD2ucoef;
//   getchar();
}
// cim5 improves on cim4, ad hoc fix for cross-derivative, 6 by 6 system
void cim5(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
          PBData &pb, GridData &grid)
{
   int r, s, t, i, m, n, sk, mid = 1, N = 2*mid;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];
   double d0[2*grid.dim];
   double ***dcoef[2*grid.dim];
   double **LU, **G, value;
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double alpha, beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;
   double jumpD1u = 0.0, ***jumpD1ucoef, jumpD1uxcoef[grid.dim], 
          jumpD1uxxcoef[grid.dim];
   double jumpD2u = 0.0, ***jumpD2ucoef, jumpD2uxcoef[grid.dim], 
          jumpD2uxxcoef[grid.dim];
   double ***D1ucoef[grid.dim], **D1uxcoef, **D1uxxcoef, ***D1uxxcoeflarge;
//   double ***D2[grid.dim][grid.dim], *D2uxcoef[grid.dim][grid.dim], 
//          *D2uxxcoef[grid.dim][grid.dim];
   double ***D2[grid.dim][3], *D2uxcoef[grid.dim][3], 
          *D2uxxcoef[grid.dim][3];
//   int sk2[grid.dim][grid.dim][4];
   int sk2[grid.dim][3][4];
   

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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
               printf("%4.16f ",getu(tindex,0,0,0.0,thesign,grid));
            }
            cout << endl;
         }
         cout << endl;
      }
   }


   double zerouxxcoef[grid.dim], zerouxcoef[grid.dim], **zerojumpuxxcoef;
   double ***zeroucoef;
   zeroucoef = matrix(N,N,N);
   D1uxcoef = matrix(grid.dim-1,grid.dim-1);
   D1uxxcoef = matrix(grid.dim-1,grid.dim-1);
   D1uxxcoeflarge = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   for (s = 0; s < grid.dim; s++) 
   {
      zerouxcoef[s] = 0.0;
      zerouxxcoef[s] = 0.0;
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      setvalarray(zeroucoef,sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;
   value = 0.0;

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   G = matrix(2*grid.dim-1,2*grid.dim-1);
   jumpD1ucoef = matrix(N,N,N);
   jumpD2ucoef = matrix(N,N,N);
   for (r = 0; r < grid.dim; r++)
   {
      D1ucoef[r] = matrix(N,N,N);
      for (s = 0; s < grid.dim; s++)
      {
         D2[r][s] = matrix(N,N,N);
         D2uxcoef[r][s] = new double[grid.dim];
         D2uxxcoef[r][s] = new double[grid.dim];
      }
   }
   for (r = 0; r < 2*grid.dim; r++)
      dcoef[r] = matrix(N,N,N);

   double ux[grid.dim];
   for (r = 0; r < grid.dim; r++)
      ux[r] = 0.0;

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
         if (getsk2(sk2[m][n],m,n,index,S,grid) && !globoldcim5)
         {
// make sure sindex is correctly equal to mid
            getD2(D2[m][n],m,n,sk2,sindex,grid);
            for (t = 0; t < grid.dim; t++)
            {
               D2uxcoef[m][n][t] = 0.0;
               D2uxxcoef[m][n][t] = 0.0;
            }
         }
         else
         {
// version with data just used this part of if statement
            getcim5D2(D2[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],m,n,index,mid,S,grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {   
               cout << "D2u " << m << " " << n << endl;
               cout << evalcoef(0.0,D2[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],index,
                                0,0,0.0,mid,thesign,grid) << " "
                    << getD2u(index,m,n,0,0,0.0,thesign,grid) << endl;
            }
         }

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
               cout << "r = " << r << " and sk = " << sk << endl;
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               for (s = 0; s < grid.dim; s++)
                  setvalarray(D1ucoef[s],sindex,0.0);
               setvalarray(jumpD1ucoef,sindex,0.0);

               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
            
            rindex[r] = index[r]+sk;
// getting interface info
            getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
            beta = 1.0-alpha;
// get derivatives
            getDu(D1uxcoef,D1uxxcoeflarge,index,r,sk,alpha,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               if (index[0] == eindex[0] && index[1] == eindex[1] && 
                   index[2] == eindex[2])
               {
                  cout << "   Du before " << m << ": ";
                  cout << evalcoef(0.0,zeroucoef,D1uxcoef[m],D1uxxcoeflarge[m],index,
                                   0,0,0.0,mid,thesign,grid) << " "
                       << getDu(index,m,r,sk,alpha,thesign,grid) << endl;
               }
            for (m = 0; m < grid.dim; m++)
               recast(D1ucoef[m],D1uxcoef[m],D1uxxcoeflarge[m],D2,D2uxcoef,D2uxxcoef,
                      mid,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = 0; n < grid.dim; n++)
                  if (globoldwrongcim5)
                     D1uxxcoef[m][n] = D1uxxcoeflarge[m][m][n];
                  else
                     D1uxxcoef[m][n] = D1uxxcoeflarge[m][n][n];
            for (m = 0; m < grid.dim; m++)
               if (index[0] == eindex[0] && index[1] == eindex[1] && 
                   index[2] == eindex[2])
               {
                  cout << "   Du " << m << ": ";
                  cout << evalcoef(0.0,D1ucoef[m],D1uxcoef[m],D1uxxcoef[m],index,0,0,
                                   0.0,mid,thesign,grid) << " "
                       << getDu(index,m,r,sk,alpha,thesign,grid) << endl;
               }
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha,grid);
            getjumpux(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,index,r,sk,alpha,
                      thesign,normal,tangent,mid,D1ucoef,D1uxcoef,D1uxxcoef,S,pb,grid);
            getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,alpha,
                       thesign,normal,mid,D1ucoef,D1uxcoef,D1uxxcoef,D2,D2uxcoef,
                       D2uxxcoef,S,pb,grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && 
                index[2] == eindex[2])
            {
               cout << "   jumpD1u: ";
               cout << evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,index,0,
                                0,0.0,mid,thesign,grid) << " "
                    << getDu(index,r,r,sk,alpha,1,grid)-getDu(index,r,r,sk,alpha,-1,grid)
                    << endl;
               cout << "   jumpD2u " << index[0] << " " << index[1] << " " << index[2] 
                    << ": ";
               cout << evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,0,
                                0,0.0,mid,thesign,grid) << " "
                    << getD2u(index,r,r,r,sk,alpha,1,grid)-
                       getD2u(index,r,r,r,sk,alpha,-1,grid) << endl;
            }
// form d0 and dcoef's rth entry 
            d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                               sk*beta*jumpD1u/grid.dx[r]+
                                               0.5*beta*beta*jumpD2u);

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
   for (j = 0; j < grid.dim; j++)
      uxx[j] = 0.0;

   gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,2*grid.dim-1);
   value = 0.0;
   for (n = 0; n < grid.dim; n++)
      value += ehere*temp[n];
//   setvalarray(b,index,evalarray(b,index)+value);
   setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid)+value);
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
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value -= ehere*temp[n];
      if (value != 0.0)
         sparse2(index,tindex,A,value,grid);
      for (j = 0; j < grid.dim; j++)
         if (evalarray(S,tindex) < 0.0)
            uxx[j] += temp[j]*getu(tindex,0,0,0.0,-1,grid);
         else
            uxx[j] += temp[j]*getu(tindex,0,0,0.0,1,grid);
      for (j = 0; j < grid.dim; j++)
         if (evalarray(S,tindex) < 0.0)
            ux[j] += temp[j+grid.dim]*getu(tindex,0,0,0.0,-1,grid);
         else
            ux[j] += temp[j+grid.dim]*getu(tindex,0,0,0.0,1,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   if (globerr < fabs(ehere*(uxx[0]+uxx[1]+uxx[2])+
                      getf(index,0,0,0.0,thesign,pb,grid)))
   {
      globerr = fabs(ehere*(uxx[0]+uxx[1]+uxx[2])+
                     getf(index,0,0,0.0,thesign,pb,grid));
      for (s = 0; s < grid.dim; s++)
         globerrvec3[s] = index[s];
   }

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
   {
      cout << ux[0] << " " << ux[1] << " " << ux[2] << endl;
      cout << uxx[0] << " " << uxx[1] << " " << uxx[2] << endl;
      // getchar();
   }
//      exit(1);
   free_matrix(zeroucoef,N,N,N);
   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(G,2*grid.dim-1,2*grid.dim-1);
   free_matrix(jumpD1ucoef,N,N,N);
   free_matrix(jumpD2ucoef,N,N,N);
   free_matrix(D1uxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1uxxcoeflarge,grid.dim-1,grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      free_matrix(D1ucoef[r],N,N,N);
      for (s = 0; s < grid.dim; s++)
      {
         free_matrix(D2[r][s],N,N,N);
         delete [] D2uxcoef[r][s];
         delete [] D2uxxcoef[r][s];
      }
   }
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(dcoef[r],N,N,N);
}
// cim345 with no dusmall
void cim345(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
            PBData &pb, GridData &grid)
{
   int r, s, t, i, m, n, sk, mid = 2, N = 2*mid;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];
   double d0[2*grid.dim];
   double ***dcoef[2*grid.dim];
   double **LU, **G, value;
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
   

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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

   double zerouxxcoef[grid.dim], zerouxcoef[grid.dim];
   double ***zeroucoef;
   zeroucoef = matrix(N,N,N);
   for (s = 0; s < grid.dim; s++) 
   {
      zerouxcoef[s] = 0.0;
      zerouxxcoef[s] = 0.0;
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      setvalarray(zeroucoef,sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;
   value = 0.0;

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   G = matrix(2*grid.dim-1,2*grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      dcoef[r] = matrix(N,N,N);

   double ux[grid.dim];
   for (r = 0; r < grid.dim; r++)
      ux[r] = 0.0;

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
//            cout << "r = " << r << " " << sk << endl;
//            cout << "getting D2u" << endl;
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
                  if (index[0] == eindex[0] && index[1] == eindex[1] && 
                      index[2] == eindex[2])
                     cout << "D2u " << m << " " << n << " " << D2u[m][n] << " " 
                          << getD2u(index,m,n,0,0,0.0,thesign,grid) << " "
                          << fabs(D2u[m][n]-getD2u(index,m,n,0,0,0.0,thesign,grid))
                          <<  endl;
               }
            
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
   for (j = 0; j < grid.dim; j++)
      uxx[j] = 0.0;

   gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,2*grid.dim-1);
   value = 0.0;
   for (n = 0; n < grid.dim; n++)
      value += ehere*temp[n];
//   setvalarray(b,index,evalarray(b,index)+value);
   setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid)+value);
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
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value -= ehere*temp[n];
      if (value != 0.0)
         sparse2(index,tindex,A,value,grid);
      for (j = 0; j < grid.dim; j++)
         if (evalarray(S,tindex) < 0.0)
            uxx[j] += temp[j]*getu(tindex,0,0,0.0,-1,grid);
         else
            uxx[j] += temp[j]*getu(tindex,0,0,0.0,1,grid);
      for (j = 0; j < grid.dim; j++)
         if (evalarray(S,tindex) < 0.0)
            ux[j] += temp[j+grid.dim]*getu(tindex,0,0,0.0,-1,grid);
         else
            ux[j] += temp[j+grid.dim]*getu(tindex,0,0,0.0,1,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   if (globerr < fabs(ehere*(uxx[0]+uxx[1]+uxx[2])+
                      getf(index,0,0,0.0,thesign,pb,grid)))
   {
      globerr = fabs(ehere*(uxx[0]+uxx[1]+uxx[2])+
                     getf(index,0,0,0.0,thesign,pb,grid));
      for (s = 0; s < grid.dim; s++)
         globerrvec3[s] = index[s];
   }
//   cout << ux[0] << " " << ux[1] << " " << ux[2] << endl;
//   cout << uxx[0] << " " << uxx[1] << " " << uxx[2] << endl;
//   getchar();

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
   {
      cout << ux[0] << " " << ux[1] << " " << ux[2] << endl;
      cout << uxx[0] << " " << uxx[1] << " " << uxx[2] << endl;
      cout << endl;
      outputneighbors(S,index,grid);
      // getchar();
   }
//      exit(1);
   free_matrix(zeroucoef,N,N,N);
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
// cim345 with Dusmall
void cim345(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, 
            int *index, int gamma[][2], double ***S, PBData &pb, GridData &grid)
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
   
   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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

            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               cout<<endl<<"dim = "<<r <<" sk = "<<sk << endl;
            }
            
            // getting interface info            
            rindex[r] = index[r]+sk;
            getinterfaceinfo(alpha[r][(sk+1)/2],tangent,normal,S,index,rindex,grid);
            beta = 1.0-alpha[r][(sk+1)/2];

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

                  #ifdef FIXBANANA 
                  //surgical fix for banana
                  if ( m==0 && n==2  && GRIDNUM==110 && SURFOPT==13 && index[0]==18 && (index[1]==54 || index[1]==56) && (index[2]==11||index[2]==99))
                  {
                      // printf("fix (%d,%d) plane at index(%d,%d,%d)\n",m,n,index[0],index[1],index[2]);
                    cout<<"fix "<<m<<" "<<n<<" at "<<index[0]<<","<<index[1]<<","<<index[2]<<endl;
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
                      }
                      int minIdx = std::min_element(err.begin(),err.end()) - err.begin();
                      
                      copyMat(D2ucoef[m][n],D2ucoefvec[minIdx],N,N,N);
                      copy(D2uxcoefvec[minIdx],D2uxcoefvec[minIdx]+3, D2uxcoef[m][n]);
                      copy(D2uxxcoefvec[minIdx],D2uxxcoefvec[minIdx]+3, D2uxxcoef[m][n]);

                      for(int i = 0; i < D2ucoefvec.size(); i++){
                        free_matrix(D2ucoefvec[i],N,N,N);
                        delete [] D2uxcoefvec[i];
                        delete [] D2uxxcoefvec[i];
                      }
                    }
                  #endif

                  if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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

// get derivatives
//            cout << "getting Du" << endl;
            getcim345Du(D1u[r][(sk+1)/2],D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                        D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,index,r,sk,
                        alpha[r][(sk+1)/2],thesign,D2u,D2ucoef,D2uxcoef,D2uxxcoef,
                        D2jumpuxxcoef,mid,grid);
            double rhs = 0.0;
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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
            getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha[r][(sk+1)/2],thesign,normal,mid,D1u[r][(sk+1)/2],
                             D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                             D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
                             D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                             jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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

   double uxval[grid.dim], uxxval[grid.dim]; // for chekcing Dusmall

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
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
   // write G matrix
   if (globwritemx){
      outfile_mx<<index[0]<<","<<index[1]<<","<<index[2]<<",";
      for (m = 0; m < 2*grid.dim; m++){
           for (n = 0; n < 2*grid.dim; n++)
              outfile_mx <<setprecision(12)<<scientific<< G[m][n] << ",";
        }
      outfile_mx<<endl;
   }

// get G [uxx .. ux] = d0 + dcoef u. 
// after inversion, we get 
// u_rr = D2u[r][r] + D2ucoef[r][r][sindex] u[sindex] 
// u_r = ux[r] +  uxcoef[r][sindex] u[sindex]
// solve for uxx and put in matrix
   gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,2*grid.dim-1);
   value = 0.0;
   for (n = 0; n < grid.dim; n++)
      value += ehere*temp[n];
//   setvalarray(b,index,evalarray(b,index)+value);
   setvalarray(b,index,getf(index,0,0,0.0,thesign,   pb,grid)+value);
   for (j = 0; j < grid.dim; j++)
      D2u[j][j] = temp[j];
   for (j = 0; j < grid.dim; j++)
      ux[j] = temp[j+grid.dim];

   // for (j = 0; j < grid.dim; j++)
   // {
   //    uxval[j] = temp[j+grid.dim];
   //    uxxval[j] = temp[j];
   // }

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

      // for (n = 0; n < grid.dim; n++)
      //    if (evalarray(S,tindex) < 0.0)
      //    {
      //       uxval[n] += temp[n+grid.dim]*getu(tindex,0,0,0.0,-1,grid);
      //       uxxval[n] += temp[n]*getu(tindex,0,0,0.0,-1,grid);
      //    }
      //    else
      //    {
      //       uxval[n] += temp[n+grid.dim]*getu(tindex,0,0,0.0,1,grid);
      //       uxxval[n] += temp[n]*getu(tindex,0,0,0.0,1,grid);
      //    }

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

// storing uint info
   double tol = 1.0e-14;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
//            cout << uxval[r] << " " << getDu(index,r,0,0,0.0,thesign,grid) << endl;
//            cout << uxxval[r] << " " << getD2u(index,r,r,0,0,0.0,thesign,grid) << endl;
//            cout << "uint = " << getu(index,0,0,0.0,thesign,grid)+
//                                 sk*alpha[r][(sk+1)/2]*grid.dx[r]*uxval[r]+
//                                 0.5*(alpha[r][(sk+1)/2]*grid.dx[r])*
//                                     (alpha[r][(sk+1)/2]*grid.dx[r])*uxxval[r]
//                 << endl;
//            cout << getu(index,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;

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

/*
            uxval[r] = ux[r];
            uxxval[r] = D2u[r][r];
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               for (s = 0; s < grid.dim; s++)
                  tindex[s] = index[s]-mid+sindex[s];
               if (evalarray(S,tindex) < 0.0)
               {
                  uxval[r] += evalarray(uxcoef[r],sindex)*getu(tindex,0,0,0.0,-1,grid);
                  uxxval[r] += evalarray(D2ucoef[r][r],sindex)*
                               getu(tindex,0,0,0.0,-1,grid);
               }
               else
               {
                  uxval[r] += evalarray(uxcoef[r],sindex)*getu(tindex,0,0,0.0,1,grid);
                  uxxval[r] += evalarray(D2ucoef[r][r],sindex)*
                               getu(tindex,0,0,0.0,1,grid);
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
            cout << uxval[r] << " " << getDu(index,r,0,0,0.0,thesign,grid) << endl;
            cout << uxxval[r] << " " << getD2u(index,r,r,0,0,0.0,thesign,grid) << endl;

            int Narray[grid.dim];
            for (s = 0; s < grid.dim; s++)
               Narray[s] = N;
            value = 0.0;
            for (SparseElt* current = Dusmall[buildsize-1].head; current != NULL; 
                 current = (*current).next)
               if ((*current).cindex >= 0)
               {
                  ind2sub(sindex,(*current).cindex,Narray,grid.dim);
                  for (s = 0; s < grid.dim; s++)
                     tindex[s] = index[s]+sindex[s]-mid;
                  if (evalarray(S,tindex) < 0.0)
                     value += (*current).val*getu(tindex,0,0,0.0,-1,grid);
                  else
                     value += (*current).val*getu(tindex,0,0,0.0,1,grid);
               }
               else
                  value += (*current).val;
            cout << "another uint = " << value << endl;
//            delta2d << index[0] << " " << index[1] << " " << index[2] << " " << r 
//                    << " " << sk << " " <<  value << endl;
//            getchar();
*/

// storing Du info
            for (t = 0; t < grid.dim; t++)
            {
//               cout << index[0] << " " << index[1] << " " << index[2] << " " 
//                    << r << " " << sk << " " << t << endl;
//               cout << "real = " << getDu(index,t,r,sk,alpha[r][(sk+1)/2],thesign,grid) 
//                    << endl;
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

//               uxval[t] = D1u[r][(sk+1)/2][t];
//               for (m = 0; m < grid.dim; m++)
//                  uxval[t] += D1uxcoef[r][(sk+1)/2][t][m]*ux[m]+
//                              D1uxxcoef[r][(sk+1)/2][t][m]*D2u[m][m];

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

//                  for (m = 0; m < grid.dim; m++)
//                     tindex[m] = index[m]-mid+sindex[m];
//                  if (evalarray(S,tindex) < 0.0)
//                  {
//                     uxval[t] += evalarray(D1ucoef[r][(sk+1)/2][t],sindex)*
//                                 getu(tindex,0,0,0.0,-1,grid);
//                     for (m = 0; m < grid.dim; m++)
//                        uxval[t] += (D1uxcoef[r][(sk+1)/2][t][m]*
//                                     evalarray(uxcoef[m],sindex)+
//                                     D1uxxcoef[r][(sk+1)/2][t][m]*
//                                     evalarray(D2ucoef[m][m],sindex))*
//                                    getu(tindex,0,0,0.0,-1,grid);
//                  }
//                  else
//                  {
//                     uxval[t] += evalarray(D1ucoef[r][(sk+1)/2][t],sindex)*
//                                 getu(tindex,0,0,0.0,1,grid);
//                     for (m = 0; m < grid.dim; m++)
//                        uxval[t] += (D1uxcoef[r][(sk+1)/2][t][m]*
//                                     evalarray(uxcoef[m],sindex)+
//                                     D1uxxcoef[r][(sk+1)/2][t][m]*
//                                     evalarray(D2ucoef[m][m],sindex))*
//                                    getu(tindex,0,0,0.0,1,grid);
//                  }

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

//               cout << "computed = " << uxval[t] << endl;
//               getchar();
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

void cim345neu(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, 
               int *index, int gamma[][2], double ***S, PBData &pb, GridData &grid)
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
            getcim345Du(D1u[r][(sk+1)/2],D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                        D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,index,r,sk,
                        alpha[r][(sk+1)/2],thesign,D2u,D2ucoef,D2uxcoef,D2uxxcoef,
                        D2jumpuxxcoef,mid,grid);
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha[r][(sk+1)/2],grid);
            getcim345jumpux(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                            jumpD1jumpuxxcoef,index,r,sk,alpha[r][(sk+1)/2],thesign,
                            normal,tangent,mid,D1ucoef[r][(sk+1)/2],
                            D1uxcoef[r][(sk+1)/2],D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,
                            S,pb,grid);
            getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha[r][(sk+1)/2],thesign,normal,mid,D1u[r][(sk+1)/2],
                             D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                             D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
                             D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                             jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);
// form d0 and dcoef's rth entry 
            d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                               sk*beta*jumpD1u/grid.dx[r]+
                                               0.5*beta*beta*jumpD2u)-
                                      sk*D1u[r][(sk+1)/2][r]/grid.dx[r];

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                           thesign*(sk*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                    0.5*beta*beta*evalarray(jumpD2ucoef,sindex))-
                           sk*evalarray(D1ucoef[r][(sk+1)/2][r],sindex)/grid.dx[r]);
   
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

   double uxval[grid.dim], uxxval[grid.dim];

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

   for (j = 0; j < grid.dim; j++)
   {
      uxval[j] = temp[j+grid.dim];
      uxxval[j] = temp[j];
   }

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

      for (n = 0; n < grid.dim; n++)
         if (evalarray(S,tindex) < 0.0)
         {
            uxval[n] += temp[n+grid.dim]*getu(tindex,0,0,0.0,-1,grid);
            uxxval[n] += temp[n]*getu(tindex,0,0,0.0,-1,grid);
         }
         else
         {
            uxval[n] += temp[n+grid.dim]*getu(tindex,0,0,0.0,1,grid);
            uxxval[n] += temp[n]*getu(tindex,0,0,0.0,1,grid);
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

// storing uint info
   double tol = 1.0e-14;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
//            cout << uxval[r] << " " << getDu(index,r,0,0,0.0,thesign,grid) << endl;
//            cout << uxxval[r] << " " << getD2u(index,r,r,0,0,0.0,thesign,grid) << endl;
//            cout << "uint = " << getu(index,0,0,0.0,thesign,grid)+
//                                 sk*alpha[r][(sk+1)/2]*grid.dx[r]*uxval[r]+
//                                 0.5*(alpha[r][(sk+1)/2]*grid.dx[r])*
//                                     (alpha[r][(sk+1)/2]*grid.dx[r])*uxxval[r]
//                 << endl;
//            cout << getu(index,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;

            Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
            Dusmall[buildsize].info[1] = r;
            Dusmall[buildsize].info[2] = sk;
            Dusmall[buildsize].info[3] = -1;
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

/*
            uxval[r] = ux[r];
            uxxval[r] = D2u[r][r];
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               for (s = 0; s < grid.dim; s++)
                  tindex[s] = index[s]-mid+sindex[s];
               if (evalarray(S,tindex) < 0.0)
               {
                  uxval[r] += evalarray(uxcoef[r],sindex)*getu(tindex,0,0,0.0,-1,grid);
                  uxxval[r] += evalarray(D2ucoef[r][r],sindex)*
                               getu(tindex,0,0,0.0,-1,grid);
               }
               else
               {
                  uxval[r] += evalarray(uxcoef[r],sindex)*getu(tindex,0,0,0.0,1,grid);
                  uxxval[r] += evalarray(D2ucoef[r][r],sindex)*
                               getu(tindex,0,0,0.0,1,grid);
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
            cout << uxval[r] << " " << getDu(index,r,0,0,0.0,thesign,grid) << endl;
            cout << uxxval[r] << " " << getD2u(index,r,r,0,0,0.0,thesign,grid) << endl;

            int Narray[grid.dim];
            for (s = 0; s < grid.dim; s++)
               Narray[s] = N;
            value = 0.0;
            for (SparseElt* current = Dusmall[buildsize-1].head; current != NULL; 
                 current = (*current).next)
               if ((*current).cindex >= 0)
               {
                  ind2sub(sindex,(*current).cindex,Narray,grid.dim);
                  for (s = 0; s < grid.dim; s++)
                     tindex[s] = index[s]+sindex[s]-mid;
                  if (evalarray(S,tindex) < 0.0)
                     value += (*current).val*getu(tindex,0,0,0.0,-1,grid);
                  else
                     value += (*current).val*getu(tindex,0,0,0.0,1,grid);
               }
               else
                  value += (*current).val;
            cout << "another uint = " << value << endl;
//            delta2d << index[0] << " " << index[1] << " " << index[2] << " " << r 
//                    << " " << sk << " " <<  value << endl;
//            getchar();
*/

// storing Du info
            for (t = 0; t < grid.dim; t++)
            {
//               cout << index[0] << " " << index[1] << " " << index[2] << " " 
//                    << r << " " << sk << " " << t << endl;
//               cout << "real = " << getDu(index,t,r,sk,alpha[r][(sk+1)/2],thesign,grid) 
//                    << endl;
               Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
               Dusmall[buildsize].info[1] = r;
               Dusmall[buildsize].info[2] = sk;
               Dusmall[buildsize].info[3] = t;

               value = D1u[r][(sk+1)/2][t];
               for (m = 0; m < grid.dim; m++)
                  value += D1uxcoef[r][(sk+1)/2][t][m]*ux[m]+
                           D1uxxcoef[r][(sk+1)/2][t][m]*D2u[m][m];
               if (fabs(value) > tol)
                  sparseorder(-1,Dusmall[buildsize].head,value);

//               uxval[t] = D1u[r][(sk+1)/2][t];
//               for (m = 0; m < grid.dim; m++)
//                  uxval[t] += D1uxcoef[r][(sk+1)/2][t][m]*ux[m]+
//                              D1uxxcoef[r][(sk+1)/2][t][m]*D2u[m][m];

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

//                  for (m = 0; m < grid.dim; m++)
//                     tindex[m] = index[m]-mid+sindex[m];
//                  if (evalarray(S,tindex) < 0.0)
//                  {
//                     uxval[t] += evalarray(D1ucoef[r][(sk+1)/2][t],sindex)*
//                                 getu(tindex,0,0,0.0,-1,grid);
//                     for (m = 0; m < grid.dim; m++)
//                        uxval[t] += (D1uxcoef[r][(sk+1)/2][t][m]*
//                                     evalarray(uxcoef[m],sindex)+
//                                     D1uxxcoef[r][(sk+1)/2][t][m]*
//                                     evalarray(D2ucoef[m][m],sindex))*
//                                    getu(tindex,0,0,0.0,-1,grid);
//                  }
//                  else
//                  {
//                     uxval[t] += evalarray(D1ucoef[r][(sk+1)/2][t],sindex)*
//                                 getu(tindex,0,0,0.0,1,grid);
//                     for (m = 0; m < grid.dim; m++)
//                        uxval[t] += (D1uxcoef[r][(sk+1)/2][t][m]*
//                                     evalarray(uxcoef[m],sindex)+
//                                     D1uxxcoef[r][(sk+1)/2][t][m]*
//                                     evalarray(D2ucoef[m][m],sindex))*
//                                    getu(tindex,0,0,0.0,1,grid);
//                  }

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

//               cout << "computed = " << uxval[t] << endl;
//               getchar();
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

double evalcoef(double u0, double ***ucoef, double *uxxcoef, int *index, double ***S, 
                GridData grid)
{
// only for mid = 1
   int i, s, tindex[grid.dim], sindex[grid.dim];
   double value = u0;

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-1+sindex[s];
      if (evalarray(S,tindex) < 0.0)
         value += evalarray(ucoef,sindex)*getu(tindex,0,0,0.0,-1,grid);
      else
         value += evalarray(ucoef,sindex)*getu(tindex,0,0,0.0,1,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;

   if (evalarray(S,index) < 0.0)
      for (s = 0; s < grid.dim; s++)
         value += uxxcoef[s]*getD2u(index,s,s,0,0,0.0,-1,grid);
   else
      for (s = 0; s < grid.dim; s++)
         value += uxxcoef[s]*getD2u(index,s,s,0,0,0.0,1,grid);

   return value;
}

double evalcoef(double u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                double **jumpuxxcoef, int *index, int rstar, int sstar, double alpha, 
                double ***S, GridData grid)
{
   int i, s, tindex[grid.dim], sindex[grid.dim];
   double value = u0;

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 5)
   {
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-2+sindex[s];
      if (evalarray(S,tindex) < 0.0)
         value += evalarray(ucoef,sindex)*getu(tindex,0,0,0.0,-1,grid);
      else
         value += evalarray(ucoef,sindex)*getu(tindex,0,0,0.0,1,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 2;

   if (evalarray(S,index) < 0.0)
      for (s = 0; s < grid.dim; s++)
         value += uxxcoef[s]*getD2u(index,s,s,0,0,0.0,-1,grid);
   else
      for (s = 0; s < grid.dim; s++)
         value += uxxcoef[s]*getD2u(index,s,s,0,0,0.0,1,grid);

   if (evalarray(S,index) < 0.0)
      for (s = 0; s < grid.dim; s++)
         value += uxcoef[s]*getDu(index,s,0,0,0.0,-1,grid);
   else
      for (s = 0; s < grid.dim; s++)
         value += uxcoef[s]*getDu(index,s,0,0,0.0,1,grid);

   for (i = 0; i < grid.dim; i++)
      for (s = i; s < grid.dim; s++)
         value += jumpuxxcoef[i][s]*(getD2u(index,i,s,rstar,sstar,alpha,1,grid)-
                                     getD2u(index,i,s,rstar,sstar,alpha,-1,grid));

   return value;
}

double evalcoef(double u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                int *index, int rstar, int sstar, double alpha, int mid, 
                double thesign, GridData grid)
{
   int i, s, N = 2*mid, tindex[grid.dim], sindex[grid.dim];
   double value = u0;

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-mid+sindex[s];
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

   for (s = 0; s < grid.dim; s++)
   {
      value += uxcoef[s]*getDu(index,s,rstar,sstar,alpha,thesign,grid);
      value += uxxcoef[s]*getD2u(index,s,s,rstar,sstar,alpha,thesign,grid);
   }

   return value;
}

double evalcoef(double u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                int *index, int rstar, int sstar, double alpha, int mid, 
                double ***S, GridData grid)
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
//uxx include cross derivative
double evalcoef(double u0, double ***ucoef, double *uxcoef, double **uxxcoef, 
                int *index, int rstar, int sstar, double alpha, int mid, 
                double thesign, GridData grid)
{
   int i, s, t, N = 2*mid, tindex[grid.dim], sindex[grid.dim];
   double value = u0;

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-mid+sindex[s];
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

   for (s = 0; s < grid.dim; s++)
   {
      value += uxcoef[s]*getDu(index,s,rstar,sstar,alpha,thesign,grid);
      for (t = 0; t < grid.dim; t++)
         value += uxxcoef[s][t]*getD2u(index,s,t,rstar,sstar,alpha,thesign,grid);
   }

   return value;
}

/*
double evalcoef(double u0, double ***ucoef, double *uxcoef, double **uxxcoef, 
                double **jumpuxxcoef, int *index, int rstar, int sstar, double alpha, 
                int mid, double thesign, GridData grid)
{
   int i, s, t, N = 2*mid, tindex[grid.dim], sindex[grid.dim];
   double value = u0;

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-mid+sindex[s];
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

   for (s = 0; s < grid.dim; s++)
   {
      value += uxcoef[s]*getDu(index,s,rstar,sstar,alpha,thesign,grid);
      for (t = 0; t < grid.dim; t++)
         value += uxxcoef[s][t]*getD2u(index,s,t,rstar,sstar,alpha,thesign,grid);
   }

   return value;
}
*/
// D1 is 3x3x3x3, D[r] is coefficient for u in 3x3x3 kernel when approximating derivative of Du in dim r
// D1 only include first order component. 
void getD1(double ****D1, int *sk, GridData grid)
{
   int i, m, s, sindex[grid.dim];

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)// set all as 0
   {
      for (m = 0; m < grid.dim; m++)
         setvalarray(D1[m],sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;
   for (m = 0; m < grid.dim; m++)
   {
      if (sk[m] == 0)// if no interface ins dim m, central difference
      {
         for (s = -1; s <= 1; s += 2)
         {
            sindex[m] = 1+s;
            setvalarray(D1[m],sindex,s/(2.0*grid.dx[m]));
         }
         sindex[m] = 1;
      }
      else //if exist interface in derection m, backward difference
      {
         for (s = 0; s <= 1; s++)
         {
            sindex[m] = 1-s*sk[m];
            setvalarray(D1[m],sindex,(1-2*s)*sk[m]/grid.dx[m]);
         }
         sindex[m] = 1;
      }
   }
}
// coefficient of cross derivative in (m,n) plane in terms of u in 3x3x3 kernel
void getD2(double ***D2[][3], int sk2[][3][4], GridData grid)
{
   int i, m, n, s, t, sindex[grid.dim];

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
            setvalarray(D2[m][n],sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
// get sk2
         for (s = 0; s < 2; s++)
         {
            sindex[m] = 1+sk2[m][n][s];
            for (t = 0; t < 2; t++)
            {
               sindex[n] = 1+sk2[m][n][2+t];
               setvalarray(D2[m][n],sindex,(2*s-1)*(2*t-1)/
                                           ((sk2[m][n][1]-sk2[m][n][0])*
                                            (sk2[m][n][3]-sk2[m][n][2])*
                                            grid.dx[m]*grid.dx[n]));
            }
            sindex[n] = 1;
         }
         sindex[m] = 1;
      }
}

void getD2(double ***D2, int m, int n, int sk2[][3][4], GridData grid)
{
   int i, s, t, sindex[grid.dim];

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      setvalarray(D2,sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;
   for (s = 0; s < 2; s++)
   {
      sindex[m] = 1+sk2[m][n][s];
      for (t = 0; t < 2; t++)
      {
         sindex[n] = 1+sk2[m][n][2+t];
         setvalarray(D2,sindex,
                     (2*s-1)*(2*t-1)/((sk2[m][n][1]-sk2[m][n][0])*
                                      (sk2[m][n][3]-sk2[m][n][2])*
                                      grid.dx[m]*grid.dx[n]));
      }
      sindex[n] = 1;
   }
   sindex[m] = 1;
}

void getD2(double ***D2, int m, int n, int sk2[][3][4], int *tindex, GridData grid)
{
   int i, s, t, N[grid.dim], sindex[grid.dim];

   for (s = 0; s < grid.dim; s++)
      N[s] = 2*tindex[s];
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
      sindex[m] = tindex[m]+sk2[m][n][s];
      for (t = 0; t < 2; t++)
      {
         sindex[n] = tindex[n]+sk2[m][n][2+t];
         setvalarray(D2,sindex,
                     (2*s-1)*(2*t-1)/((sk2[m][n][1]-sk2[m][n][0])*
                                      (sk2[m][n][3]-sk2[m][n][2])*
                                      grid.dx[m]*grid.dx[n]));
      }
      sindex[n] = tindex[n];
   }
   sindex[m] = tindex[m];
}
//used in CIM4, cim12.pdf page 36, D2[5x5x5] is cross-derivative in m,n plane in terms of u-value
void getD2(double ***D2, double &jumpuxxcoef, int m, int n, int rstar, int sstar, 
           double thesign, int sk2[][3][4], GridData grid)
{
   int i, s, t, sindex[grid.dim], tindex[grid.dim];

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 5)
   {
      setvalarray(D2,sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }

   for (s = 0; s < grid.dim; s++)
      tindex[s] = 2;
   tindex[rstar] = 2+sstar; //tindex is point across interface
   for (s = 0; s < grid.dim; s++)
      sindex[s] = tindex[s];
   for (s = 0; s < 2; s++)
   {
      sindex[m] = tindex[m]+sk2[m][n][s];
      for (t = 0; t < 2; t++)
      {
         sindex[n] = tindex[n]+sk2[m][n][2+t];
         setvalarray(D2,sindex,
                     (2*s-1)*(2*t-1)/((sk2[m][n][1]-sk2[m][n][0])*
                                      (sk2[m][n][3]-sk2[m][n][2])*
                                      grid.dx[m]*grid.dx[n]));
      }
      sindex[n] = tindex[n];
   }
   sindex[m] = tindex[m];

   jumpuxxcoef = thesign;
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
         setvalarray(D2,sindex,
                     (2*s-1)*(2*t-1)/((sk2[1]-sk2[0])*(sk2[3]-sk2[2])*
                                      grid.dx[m]*grid.dx[n]));
      }
      sindex[n] = tindex[n];
   }
   sindex[m] = tindex[m];
}
// approximate D2u in (m,n) plane at index
// in terms of constant u0, NxNxN u coeff, 3x1 uxcoeff, 3x1 uxx coef, 1x1 jumpuxx coef 
// In the first pass, rstar=0, sstar=0, return perm=1 if can be approx with same side, otherwise perm=0
// if the first pass return 0, in the second pass, provide extra info on interface location, see if use uxy on the other side (cim4)
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
      else if (getcim5D2(u,uxcoef,uxxcoef,m,n,index,mid,S,grid))
      {
         u0 = 0.0;
         jumpuxxcoef = 0.0;
         perm = 1;
      }
      else 
      {
         int r, s;
         int **anosk2 = imatrix(1,3);//2x4 matrix
         char yes[2];

         if (globdirectD2) // approximate D2 from nbr of same side, prefer out of plane
         {
            for (t = 0; t < grid.dim; t++)
               rindex[t] = index[t];
            for (r = 0; r < grid.dim && (r == m || r == n); r++); // find r no equal to m or n
            for (s = -1; s <= 1; s += 2)
            {
               rindex[r] = index[r]+s; //rindex is nbr in one of m,n-plane
               if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) != 1)//if change sign
                  yes[(s+1)/2] = getsk2(anosk2[(s+1)/2],m,n,rindex,S,grid); 
                  //yes[0] = 1 if s=-1, rindex has approx of mixed derivative
               else
                  yes[(s+1)/2] = 0;
            }
         }
         if (globdirectD2 && (yes[0] || yes[1])) // if we could found s = +/-1
         {  
            if ((yes[0] && !(yes[1])) || (yes[0] && yes[1] && 
                                          abs(anosk2[0][0]-anosk2[0][1])+
                                          abs(anosk2[0][2]-anosk2[0][3]) >
                                          abs(anosk2[1][0]-anosk2[1][1])+
                                          abs(anosk2[1][2]-anosk2[1][3])))
              // s=-1 has sk2 and s=1 no sk2, or 
            // s=-1 has sk2, s=1 has sk2, s=-1 has larger stencil,i.e. prefer central differencing
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
         else
            perm = 0;
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
// in the first pass, rstart = sstar = 0, return perm=2 if there is usual mix Du, then second pass is skipped
// in the second pass, location for interface is provided, 
// if perm=1 in the first pass, that is, no cim5 cim3 available, then use cim4 in rstar dim and sstar direction
// with globdist = 1, globdistvar = 1
// regsk2 with central diff > cim5 > regsk2 without central diff > sk2 on neighboring point(prefer the side with central diff) > sk2 on the other side
// with globdist = 1, globdistvar = 0, this is the same as globdist = 0, getcim345D2
// regsk2 (with or without central diff) > cim 5 > sk2 on neighboring point(prefer the side with central diff) > sk2 on the other side
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
      if (perm != 1) //perm = 0
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
// abs(sk2reg[1]-sk2reg[0])+abs(sk2reg[3]-sk2reg[2]) >= 3 means one direction of cross derivative use central difference
// (useuxD2 && (fabs(uxxcoef[m]) > tol || fabs(uxxcoef[n]) > tol) means cim5D2 use uxxcoef in m n plane
// in either case below , use sk2reg
//case 1: (perm != 1 && !globdistvar && regD2 && ( abs(sk2reg[1]-sk2reg[0])+abs(sk2reg[3]-sk2reg[2]) >= 3 ||(useuxD2 && (fabs(uxxcoef[m]) > tol || fabs(uxxcoef[n]) > tol))))
//case 2: (perm != 1 &&  globdistvar && regD2 &&  abs(sk2reg[1]-sk2reg[0])+abs(sk2reg[3]-sk2reg[2]) >= 3)
// if globdistvar = 0, we choose sk2reg even when sk2reg is not central diff,
// if globdistvar = 1, we choose sk2reg only when sk2reg is central diff     

      if ((perm != 1 && !globdistvar && regD2 && 
           (abs(sk2reg[1]-sk2reg[0])+abs(sk2reg[3]-sk2reg[2]) >= 3 ||
            (useuxD2 && (fabs(uxxcoef[m]) > tol || fabs(uxxcoef[n]) > tol)))) ||
          (perm != 1 && globdistvar && regD2 && 
           abs(sk2reg[1]-sk2reg[0])+abs(sk2reg[3]-sk2reg[2]) >= 3))
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
      else if (perm != 1 && (nearD2[0] || nearD2[1])) //if (regD2 and getcim5D2 failed), use nearby points, prefer central diff
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
      else if (rstar >= 0 && rstar < grid.dim && sstar != 0 && otherD2) //everything above fail, use the other side
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
// calculate [ux] in dim rstar,
// u0 is constant part, ucoef is coefficient for u, uxxcoef is coefficient for uxx
void getjumpux(double& u0, double ***ucoef, double *uxxcoef, int *index, int rstar, 
               int *sk, double alpha, int thesign, double *normal, double *tangent, 
               double ****D1ucoef, double **D1uxxcoef, double ***S, PBData &pb, 
               GridData &grid)
{
   double ehere, ethere, sigma, Dtau[grid.dim], value;
   int i, s, n, sindex[grid.dim];

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
   getsigma(sigma,index,rstar,sk[rstar],alpha,normal,pb,grid);
   getDtau(Dtau,index,rstar,sk[rstar],alpha,grid);

   u0 = sigma/ethere*normal[rstar]+getdotprod(Dtau,tangent,grid.dim)*tangent[rstar];//

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += evalarray(D1ucoef[n],sindex)*normal[n]; // dot(grad(u-),n), this part work on u coeff of u-
      setvalarray(ucoef,sindex,thesign*(ethere-ehere)/ethere*value*normal[rstar]);//[epsilon] dot(grad(u-),n) n[r] / epsp

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;

   for (s = 0; s < grid.dim; s++)
   {
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s]*normal[n]; // dot(grad(u-),n), this part work on uxx coeff of u-
      uxxcoef[s] = thesign*(ethere-ehere)/ethere*value*normal[rstar]; //[epsilon] dot(grad(u-),n) n[r] / epsp
   }
}
// in cim4, cim12.pdf page 27
void getjumpux(double& u0, double *uxcoef, double **uxxcoef, int *index, int rstar, 
               int sk, double alpha, int thesign, double *normal, double *tangent, 
               double **D1uxcoef, double ***D1uxxcoef, double ***S, PBData &pb, 
               GridData &grid)
{
   double ehere, ethere, sigma, Dtau[grid.dim], value;
   int s, t, m;

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

   u0 = sigma/ethere*normal[rstar]+getdotprod(Dtau,tangent,grid.dim)*tangent[rstar];

   for (s = 0; s < grid.dim; s++)
   {
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         value += D1uxcoef[m][s]*normal[m];//dot(grad(u-),normal)
      uxcoef[s] = thesign*(ethere-ehere)/ethere*value*normal[rstar];
      for (t = s; t < grid.dim; t++)
      {
         value = 0.0;
         for (m = 0; m < grid.dim; m++)
            value += D1uxxcoef[m][s][t]*normal[m];
         uxxcoef[t][s] = 0.0;
         uxxcoef[s][t] = thesign*(ethere-ehere)/ethere*value*normal[rstar];
      }
   }
}

void getjumpux(double& u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
               int *index, int rstar, int sk, double alpha, int thesign, 
               double *normal, double *tangent, int mid, double ****D1ucoef, 
               double **D1uxcoef, double **D1uxxcoef, double ***S, PBData &pb, 
               GridData &grid)
{
   double ehere, ethere, sigma, Dtau[grid.dim], value;
   int i, s, t, m, N = 2*mid, sindex[grid.dim];

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

   u0 = sigma/ethere*normal[rstar]+getdotprod(Dtau,tangent,grid.dim)*tangent[rstar];

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         value += evalarray(D1ucoef[m],sindex)*normal[m];
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
      for (m = 0; m < grid.dim; m++)
         value += D1uxcoef[m][s]*normal[m];
      uxcoef[s] = thesign*(ethere-ehere)/ethere*value*normal[rstar];
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         value += D1uxxcoef[m][s]*normal[m];
      uxxcoef[s] = thesign*(ethere-ehere)/ethere*value*normal[rstar];
   }
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
// used in cim4, [Du] may use [DDu], collapse into u, ux, uxx
void recast(double &u0, double ***ucoef, double *uxcoef, double **uxxcoef, 
            double **jumpD2u, double *****jumpD2ucoef, double ***jumpD2uxcoef, 
            double ***jumpD2uxxcoef, double ***D2[][3], double D2jumpuxxcoef[][3],
            GridData &grid)
{
   int i, m, n, s, sindex[grid.dim];

   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 5) 
         {//uxxcoef[m][n] DDu_mn, DDu_mn = ... D2[m][n][sindex] u(sindex)
            setvalarray(ucoef,sindex,evalarray(ucoef,sindex)+
                                     uxxcoef[m][n]*evalarray(D2[m][n],sindex)); 
            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 2;
          // if some DDu_mn need info from the other side
         if (D2jumpuxxcoef[m][n] != 0.0 && uxxcoef[m][n] != 0.0)
         {
            u0 += uxxcoef[m][n]*D2jumpuxxcoef[m][n]*jumpD2u[m][n];

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] < 5) 
            {//[Du] = ... uxxcoef[m][n] DDu_mn, DDu_mn = ... D2jumpuxxcoef[m][n] [DDu], [DDu] = ... jumpD2ucoef[m][n][sindex] u(sindex)
               setvalarray(ucoef,sindex,
                           evalarray(ucoef,sindex)+
                           uxxcoef[m][n]*D2jumpuxxcoef[m][n]*
                           evalarray(jumpD2ucoef[m][n],sindex));
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 2;
            //[Du] = ... uxxcoef[m][n] DDu_mn, DDu_mn = ... D2jumpuxxcoef[m][n] [DDu], [DDu] = ... jumpD2uxcoef[m][n][s] ux(s)
            for (s = 0; s < grid.dim; s++)
               uxcoef[s] += uxxcoef[m][n]*D2jumpuxxcoef[m][n]*jumpD2uxcoef[m][n][s];
             //[Du] = ... uxxcoef[m][n] DDu_mn, DDu_mn = ... D2jumpuxxcoef[m][n] [DDu], [DDu] = ... jumpD2uxxcoef[m][n][s] uxx(s)
            for (s = 0; s < grid.dim; s++)
               uxxcoef[s][s] += uxxcoef[m][n]*D2jumpuxxcoef[m][n]*jumpD2uxxcoef[m][n][s];
         }
         uxxcoef[m][n] = 0.0;
         uxxcoef[n][m] = 0.0;
      }
}

void recast(double ***ucoef, double *uxcoef, double **uxxcoef, double ***D2[][3], 
            double *D2uxcoef[][3], double *D2uxxcoef[][3], int mid, GridData &grid)
{
   int i, m, n, s, N = 2*mid, sindex[grid.dim];

   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] <= N)
         {
            setvalarray(ucoef,sindex,evalarray(ucoef,sindex)+
                                     uxxcoef[m][n]*evalarray(D2[m][n],sindex));
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
            uxcoef[s] += uxxcoef[m][n]*D2uxcoef[m][n][s];
         for (s = 0; s < grid.dim; s++)
            uxxcoef[s][s] += uxxcoef[m][n]*D2uxxcoef[m][n][s];

         uxxcoef[m][n] = 0.0;
         uxxcoef[n][m] = 0.0;
      }
}

void getjumpuxx(double& u0, double ***ucoef, double *uxxcoef, int *index, int rstar, 
                int *sk, double alpha, int thesign, double *normal, double ****D1ucoef,
                double **D1uxxcoef, double ****D1, double*** D2[][3], double ***S, 
                PBData &pb, GridData &grid)
{
   int i, s, n, m;
   int rindex[grid.dim], sindex[grid.dim];
   double b0[2*grid.dim], ***bcoef[2*grid.dim], c0[2*grid.dim];
   double **LU, **B, **Dn; 
   double Dndott[grid.dim], Dndots[grid.dim], dotDndot[grid.dim];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double tangent1[grid.dim], tangent2[grid.dim];
   double ethere, ehere, jumpfe;
   double value, temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim];
   double Dtau[grid.dim], **D2tau, dotD2taudot[grid.dim], Dtaudot[grid.dim];
   char theorder = 2;

   double realjumpuxx[grid.dim], realux[grid.dim], dotrealuxxdot[2];

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   B = matrix(2*grid.dim-1,2*grid.dim-1);
   Dn = matrix(grid.dim-1,grid.dim-1);
   D2tau = matrix(grid.dim-1,grid.dim-1);
   for (m = 0; m < 2*grid.dim; m++)
      bcoef[m] = matrix(2,2,2);

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

//   getinterfaceinfo(alpha,tangent1,tangent2,normal,sigma,Dn,Dsigma,jumpfe,S,index,
//                    rindex,pb,grid);
   getinterfaceinfo(tangent1,tangent2,sigma,Dn,Dsigma,jumpfe,index,rstar,sk[rstar],
                    alpha,S,pb,grid);
// form matrix
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

// form Dn dot product with various vectors
   getDtau(Dtau,index,rstar,sk[rstar],alpha,grid);
   getD2tau(D2tau,index,rstar,sk[rstar],alpha,grid);
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
// form right hand side vector: b0 without u's 
   for (n = 0; n < grid.dim; n++)
      b0[n] = (sigma/ethere-Dtaudot[0])*dotDndot[n]+dotD2taudot[n];
   b0[grid.dim] = getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                  Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2];
   b0[grid.dim+1] = getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                    Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1];
   b0[grid.dim+2] = -jumpfe;
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      for (n = 0; n < 2*grid.dim; n++)
         setvalarray(bcoef[n],sindex,0.0);
      for (n = 0; n < grid.dim; n++)
         for (m = 0; m < grid.dim; m++)
         {
            if (theorder == 2)
               setvalarray(bcoef[n],sindex,evalarray(bcoef[n],sindex)+
                                           thesign*(ethere-ehere)/ethere*
                                           dotDndot[n]*normal[m]*
                                           evalarray(D1ucoef[m],sindex));
            else
               setvalarray(bcoef[n],sindex,evalarray(bcoef[n],sindex)+
                                           thesign*(ethere-ehere)/ethere*
                                            dotDndot[n]*normal[m]*
                                            evalarray(D1[m],sindex));
         }
      for (m = 0; m < grid.dim; m++)
      {
         if (theorder == 2)
         {
            setvalarray(bcoef[grid.dim],sindex,
                        evalarray(bcoef[grid.dim],sindex)+
                        thesign*(ethere-ehere)/ethere*Dndott[m]*
                        evalarray(D1ucoef[m],sindex));
            setvalarray(bcoef[grid.dim+1],sindex,
                        evalarray(bcoef[grid.dim+1],sindex)+
                        thesign*(ethere-ehere)/ethere*Dndots[m]*
                        evalarray(D1ucoef[m],sindex));
         }
         else
         {
            setvalarray(bcoef[grid.dim],sindex,
                        evalarray(bcoef[grid.dim],sindex)+
                        thesign*(ethere-ehere)/ethere*Dndott[m]*
                        evalarray(D1[m],sindex));
            setvalarray(bcoef[grid.dim+1],sindex,
                        evalarray(bcoef[grid.dim+1],sindex)+
                        thesign*(ethere-ehere)/ethere*Dndots[m]*
                        evalarray(D1[m],sindex));
         }
         for (n = 0; n < grid.dim; n++)
            if (n != m)
            {
               setvalarray(bcoef[grid.dim],sindex,
                           evalarray(bcoef[grid.dim],sindex)+
                           thesign*(ethere-ehere)/ethere*
                           normal[m]*tangent1[n]*
                           evalarray(D2[min(m,n)][max(m,n)],sindex));
               setvalarray(bcoef[grid.dim+1],sindex,
                           evalarray(bcoef[grid.dim+1],sindex)+
                           thesign*(ethere-ehere)/ethere*
                           normal[m]*tangent2[n]*
                           evalarray(D2[min(m,n)][max(m,n)],sindex));
            }
      }
      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;
   gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);
   for (m = 0; m < grid.dim; m++)
   {
      if (theorder == 2)
      {
         for (n = 0; n < grid.dim; n++)
            c0[n] = 0.0;
         c0[grid.dim] = thesign*(ethere-ehere)/ethere*normal[m]*tangent1[m];
         c0[grid.dim+1] = thesign*(ethere-ehere)/ethere*normal[m]*tangent2[m];
         for (s = 0; s < grid.dim; s++)
         {
            for (n = 0; n < grid.dim; n++)
               c0[n] += thesign*(ethere-ehere)/ethere*dotDndot[n]*
                        normal[s]*D1uxxcoef[s][m];
            c0[grid.dim] += thesign*(ethere-ehere)/ethere*Dndott[s]*D1uxxcoef[s][m];
            c0[grid.dim+1] += thesign*(ethere-ehere)/ethere*Dndots[s]*D1uxxcoef[s][m];
         }
      }
      else
      {
         for (n = 0; n < grid.dim; n++)
            c0[n] = 0.0;
         c0[grid.dim] = thesign*(ethere-ehere)/ethere*normal[m]*tangent1[m];
         c0[grid.dim+1] = thesign*(ethere-ehere)/ethere*normal[m]*tangent2[m];
      }
      c0[grid.dim+2] = 0.0;
      for (n = 0; n < 2*grid.dim; n++)
         temp[n] = c0[n];
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      uxxcoef[m] = temp[rstar];
   }
   forwardbacksub0(temp,b0,LU,PLR,PLC,2*grid.dim-1);
   u0 = temp[rstar];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      for (n = 0; n < 2*grid.dim; n++)
         temp[n] = evalarray(bcoef[n],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      setvalarray(ucoef,sindex,temp[rstar]);
      
      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;


   cout << evalcoef(u0,ucoef,uxxcoef,index,S,grid) << " " 
        << getD2u(index,rstar,rstar,rstar,sk[rstar],alpha,1,grid)-
           getD2u(index,rstar,rstar,rstar,sk[rstar],alpha,-1,grid) << endl;

   for (n = 0; n < grid.dim; n++)
      realjumpuxx[n] = getD2u(index,n,n,rstar,sk[rstar],alpha,1,grid)-
                       getD2u(index,n,n,rstar,sk[rstar],alpha,-1,grid);
   realjumpuxx[grid.dim] = getD2u(index,0,1,rstar,sk[rstar],alpha,1,grid)-
                           getD2u(index,0,1,rstar,sk[rstar],alpha,-1,grid);
   realjumpuxx[grid.dim+1] = getD2u(index,1,2,rstar,sk[rstar],alpha,1,grid)-
                             getD2u(index,1,2,rstar,sk[rstar],alpha,-1,grid);
   realjumpuxx[grid.dim+2] = getD2u(index,0,2,rstar,sk[rstar],alpha,1,grid)-
                             getD2u(index,0,2,rstar,sk[rstar],alpha,-1,grid);
   for (n = 0; n < grid.dim; n++)
      realux[n] = getDu(index,n,rstar,sk[rstar],alpha,thesign,grid);
   dotrealuxxdot[0] = 0.0;
   dotrealuxxdot[1] = 0.0;
   for (n = 0; n < grid.dim; n++)
      for (m = 0; m < grid.dim; m++)
      {
         dotrealuxxdot[0] += normal[n]*
                             getD2u(index,n,m,rstar,sk[rstar],alpha,thesign,grid)*
                             tangent1[m];
         dotrealuxxdot[1] += normal[n]*
                             getD2u(index,n,m,rstar,sk[rstar],alpha,thesign,grid)*
                             tangent2[m];
      }
   cout << index[0] << " " << index[1] << " " << index[2] << endl;
   cout << "   checking " << getdotprod(B[0],realjumpuxx,2*grid.dim) << " "
        << (sigma/ethere+
            thesign*(ethere-ehere)/ethere*getdotprod(realux,normal,grid.dim)-
            Dtaudot[0])*dotDndot[0]+dotD2taudot[0] << endl;
   cout << "   checking " << getdotprod(B[1],realjumpuxx,2*grid.dim) << " "
        << (sigma/ethere+
            thesign*(ethere-ehere)/ethere*getdotprod(realux,normal,grid.dim)-
            Dtaudot[0])*dotDndot[1]+dotD2taudot[1] << endl;
   cout << "   checking " << getdotprod(B[2],realjumpuxx,2*grid.dim) << " "
        << (sigma/ethere+
            thesign*(ethere-ehere)/ethere*getdotprod(realux,normal,grid.dim)-
            Dtaudot[0])*dotDndot[2]+dotD2taudot[2] << endl;
   cout << "   checking " << getdotprod(B[3],realjumpuxx,2*grid.dim) << " "
        << (getdotprod(Dsigma,tangent1,grid.dim)+
            thesign*(ethere-ehere)*(dotrealuxxdot[0]+
                                    getdotprod(realux,Dndott,grid.dim)))/ethere-
            Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2] << endl;
   cout << "   checking " << getdotprod(B[4],realjumpuxx,2*grid.dim) << " "
        << (getdotprod(Dsigma,tangent2,grid.dim)+
            thesign*(ethere-ehere)*(dotrealuxxdot[1]+
                                    getdotprod(realux,Dndots,grid.dim)))/ethere-
            Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1] << endl;
   cout << "   checking " << getdotprod(B[5],realjumpuxx,2*grid.dim) << " "
        << -jumpfe << endl;

   getchar();


   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
   for (m = 0; m < 2*grid.dim; m++)
      free_matrix(bcoef[m],2,2,2);
}
// first form the system G [[uxx];[uyy];[uzz];[uxy];[uyz];[uxz]] = b0 + bcoeff[3x3x3] u + bxxcoef[1x3] [uxx;uyy;uzz]
// used in cim3again, uxx along dim rstar = u0 + ucoef[3x3x3] (u)ij in + uxxcoef[3] (u_xx)ij (u_yy)ij (u_zz)ij
void getjumpuxx(double& u0, double ***ucoef, double *uxxcoef, int *index, int rstar, 
                int *sk, double alpha, int thesign, double *normal, double ****D1ucoef,
                double **D1uxxcoef, double*** D2[][3], double ***S, PBData &pb, 
                GridData &grid)
{
   int i, r, s, n, m;
   int rindex[grid.dim], sindex[grid.dim];
   double b0[2*grid.dim], bcoef[2*grid.dim], bxxcoef[2*grid.dim];
   double **LU, **B, **Dn; 
   double Dndott[grid.dim], Dndots[grid.dim], dotDndot[grid.dim];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double tangent1[grid.dim], tangent2[grid.dim];
   double ethere, ehere, jumpfe;
   double value, temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim];
   double Dtau[grid.dim], **D2tau, dotD2taudot[grid.dim], Dtaudot[grid.dim];

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   B = matrix(2*grid.dim-1,2*grid.dim-1);
   Dn = matrix(grid.dim-1,grid.dim-1);
   D2tau = matrix(grid.dim-1,grid.dim-1);

//   double bb0[2*grid.dim];
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

   getinterfaceinfo(tangent1,tangent2,sigma,Dn,Dsigma,jumpfe,index,rstar,sk[rstar],
                    alpha,S,pb,grid);
// form matrix
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
//   for (m = 0; m < 2*grid.dim; m++)
//   {
//      for (n = 0; n < 2*grid.dim; n++)
//         cout << B[m][n] << " ";
//      cout << endl;
//   }
   gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);
            
// form Dn dot product with various vectors
   getDtau(Dtau,index,rstar,sk[rstar],alpha,grid);
   getD2tau(D2tau,index,rstar,sk[rstar],alpha,grid);
   for (n = 0; n < grid.dim; n++)
   {
      Dndott[n] = 0.0;
      Dndots[n] = 0.0;
      dotDndot[n] = 0.0; //s \grad n s
      dotD2taudot[n] = 0.0; //s \grad^2 \tau s
      Dtaudot[n] = 0.0;
   }
   for (n = 0; n < grid.dim; n++)
   {
      for (m = 0; m < grid.dim; m++)
      {
         Dndott[n] += Dn[n][m]*tangent1[m];
         Dndots[n] += Dn[n][m]*tangent2[m];
         dotD2taudot[0] += tangent1[n]*D2tau[n][m]*tangent1[m];//t1 \grad^2(t) t1
         dotD2taudot[1] += tangent2[n]*D2tau[n][m]*tangent2[m];//t2 \grad^2(t) t2
         dotD2taudot[2] += tangent1[n]*D2tau[n][m]*tangent2[m];//t1 \grad^2(t) t2
      }
      dotDndot[0] += tangent1[n]*Dndott[n]; //t1 grad(n) t1
      dotDndot[1] += tangent2[n]*Dndots[n];//t2 grad(n) t2
      dotDndot[2] += tangent1[n]*Dndots[n];//t1 grad(n) t2 = t2 grad(n) t1
      Dtaudot[0] += Dtau[n]*normal[n];// grad(tau) n
      Dtaudot[1] += Dtau[n]*tangent1[n];// grad(tau) t1
      Dtaudot[2] += Dtau[n]*tangent2[n];// grad(tau) t2
   }
// form right hand side vector: b0 without u's 
//   double ****bbcoef = matrix(2*grid.dim-1,2,2,2);
//   double tempval[2*grid.dim];
//   int tindex[grid.dim];
   for (m = 0; m < grid.dim; m++)
      b0[m] = (sigma/ethere-Dtaudot[0])*dotDndot[m]+dotD2taudot[m]; // why Dtaudot[0]
   b0[grid.dim] = getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                  Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2];
   b0[grid.dim+1] = getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                    Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1];
   b0[grid.dim+2] = -jumpfe;
//   for (m = 0; m < 2*grid.dim; m++)
//      tempval[m] = b0[m];
//   for (m = 0; m < 2*grid.dim; m++)
//      cout << tempval[m] << " ";
//   cout << endl;

   forwardbacksub0(temp,b0,LU,PLR,PLC,2*grid.dim-1);
   u0 = temp[rstar];
//   for (m = 0; m < 2*grid.dim; m++)
//      bb0[m] = temp[m];
// get coefs 
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         value += evalarray(D1ucoef[m],sindex)*normal[m];
      for (m = 0; m < grid.dim; m++)
         bcoef[m] = thesign*(ethere-ehere)/ethere*value*dotDndot[m]; //from [epsilon]/epsp grad(u-) n s_i grad(n) s_j
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         value += evalarray(D1ucoef[m],sindex)*Dndott[m];
      bcoef[grid.dim] = thesign*(ethere-ehere)/ethere*value; //row 4, from [epsilon]/epsp grad(u-) grad(n) s1, row num start from 1
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         value += evalarray(D1ucoef[m],sindex)*Dndots[m];
      bcoef[grid.dim+1] = thesign*(ethere-ehere)/ethere*value; //row 5, from [epsilon]/epsp grad(u-) grad(n) s2
      bcoef[grid.dim+2] = 0.0;

      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
            value += (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                     evalarray(D2[m][n],sindex);
      bcoef[grid.dim] += thesign*(ethere-ehere)/ethere*value; //row4 from [epsilon]/epsp  n grad^2(u-) s1
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
            value += (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                     evalarray(D2[m][n],sindex);
      bcoef[grid.dim+1] += thesign*(ethere-ehere)/ethere*value; //row 5 from [epsilon]/epsp  n grad^2(u-) s2
//      for (m = 0; m < grid.dim; m++)
//         tindex[m] = index[m]-1+sindex[m];
//      for (m = 0; m < 2*grid.dim; m++)
//         tempval[m] += bcoef[m]*getu(tindex,0,0,0.0,thesign,grid);

      forwardbacksub0(temp,bcoef,LU,PLR,PLC,2*grid.dim-1);
      setvalarray(ucoef,sindex,temp[rstar]);
//      for (m = 0; m < 2*grid.dim; m++)
//         setvalarray(bbcoef[m],sindex,temp[m]);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;

//   double **bbxxcoef = matrix(2*grid.dim-1,grid.dim-1);
   for (s = 0; s < grid.dim; s++)
   {
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         value += D1uxxcoef[m][s]*normal[m];
      for (m = 0; m < grid.dim; m++)
         bxxcoef[m] = thesign*(ethere-ehere)/ethere*value*dotDndot[m];
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         value += D1uxxcoef[m][s]*Dndott[m];
      bxxcoef[grid.dim] = thesign*(ethere-ehere)/ethere*(value+normal[s]*tangent1[s]);
      value = 0.0;
      for (m = 0; m < grid.dim; m++)
         value += D1uxxcoef[m][s]*Dndots[m];
      bxxcoef[grid.dim+1] = thesign*(ethere-ehere)/ethere*(value+normal[s]*tangent2[s]);
      bxxcoef[grid.dim+2] = 0.0;
//      for (m = 0; m < 2*grid.dim; m++)
//         tempval[m] += bxxcoef[m]*getD2u(index,s,s,0,0,0.0,thesign,grid);

      forwardbacksub0(temp,bxxcoef,LU,PLR,PLC,2*grid.dim-1);
      uxxcoef[s] = temp[rstar];
//      for (m = 0; m < 2*grid.dim; m++)
//         bbxxcoef[m][s] = temp[m];
   }
//   cout << "b val ";
//   for (m = 0; m < 2*grid.dim; m++)
//      cout << tempval[m] << " ";
//   cout << endl;

//   cout << evalcoef(bb0[0],bbcoef[0],bbxxcoef[0],index,S,grid) << " "
//        << evalcoef(bb0[1],bbcoef[1],bbxxcoef[1],index,S,grid) << " "
//        << evalcoef(bb0[2],bbcoef[2],bbxxcoef[2],index,S,grid) << endl;
//   getchar();

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
}

void getjumpuxx(double& u0, double ***ucoef, double *uxxcoef, int *index, int rstar, 
                int *sk, double alpha, int thesign, double *normal, double ****D1,
                double*** D2[][3], double ***S, PBData &pb, GridData &grid)
{
   int m, n;
   double **zerouxxcoef = matrix(grid.dim-1,grid.dim-1);

   for (m = 0; m < grid.dim; m++)
      for (n = 0; n < grid.dim; n++)
         zerouxxcoef[m][n] = 0.0;

   getjumpuxx(u0,ucoef,uxxcoef,index,rstar,sk,alpha,thesign,normal,D1,zerouxxcoef,D2,
              S,pb,grid);

   free_matrix(zerouxxcoef,grid.dim-1,grid.dim-1);
}
//used in cim4, uxx along dim rstar = u0[m][n] + ucoef[m][n][5x5x5] (u)ij in + uxxcoef[m][n][3] (u_xx)ij (u_yy)ij (u_zz)ij 
// [m][n] correspond to the row of matrix associated with [DDu_mn]
void getjumpuxx(double **u0, double *****ucoef, double ***uxcoef, double ***uxxcoef, 
                int *index, int rstar, int sk, double alpha, int thesign, 
                double *normal, double **D1uxcoef, double ***D1uxxcoef, 
                double*** D2[][3], double D2jumpuxxcoef[][3], double ***S, 
                PBData &pb, GridData &grid)
{
   int i, r, s, n, m;
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

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   B = matrix(2*grid.dim-1,2*grid.dim-1);
   Dn = matrix(grid.dim-1,grid.dim-1);
   D2tau = matrix(grid.dim-1,grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      bcoef[r] = matrix(4,4,4);

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
   }//cim12.pdf p28 p36, Du- might involve [u_xy], bring to LHS
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
         if (D2jumpuxxcoef[m][n] != 0.0)
         {
            value = 0.0;
            for (s = 0; s < grid.dim; s++)
               value += D1uxxcoef[s][m][n]*normal[s];
            for (s = 0; s < grid.dim; s++)
               B[s][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value*
                                      D2jumpuxxcoef[m][n]*dotDndot[s];
            value = 0.0;
            for (s = 0; s < grid.dim; s++)
               value += D1uxxcoef[s][m][n]*Dndott[s];
            B[grid.dim][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*
                                          (value+normal[m]*tangent1[n]+
                                                 normal[n]*tangent1[m])*
                                          D2jumpuxxcoef[m][n];
            value = 0.0;
            for (s = 0; s < grid.dim; s++)
               value += D1uxxcoef[s][m][n]*Dndots[s];
            B[grid.dim+1][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*
                                            (value+normal[m]*tangent2[n]+
                                                   normal[n]*tangent2[m])*
                                            D2jumpuxxcoef[m][n];
         }
//   double tempval[2*grid.dim], ux[grid.dim];
//   int tindex[grid.dim];
// form right hand side vector: b0 without u's 
   for (n = 0; n < grid.dim; n++)
      b0[n] = (sigma/ethere-Dtaudot[0])*dotDndot[n]+dotD2taudot[n];
   b0[grid.dim] = getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                  Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2];
   b0[grid.dim+1] = getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                    Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1];
   b0[grid.dim+2] = -jumpfe;
//   for (m = 0; m < 2*grid.dim; m++)
//      tempval[m] = b0[m];
//   ux[rstar] = sk*(getu(index,0,0,0.0,thesign,grid)-
//                   getu(index,rstar,-sk,1.0,thesign,grid))/grid.dx[rstar]+
//               sk*grid.dx[rstar]/2.0*getD2u(index,rstar,rstar,0,0,0.0,thesign,grid);
//   for (m = 0; m < grid.dim; m++)
//      if (m != rstar)
//         ux[m] = (getu(index,m,1,1.0,thesign,grid)-
//                  getu(index,m,-1,1.0,thesign,grid))/(2.0*grid.dx[rstar]);
//   cout << ux[0] << " " << ux[1] << " " << ux[2] << endl;
// get coefs
   for (s = 0; s < grid.dim; s++)
   {
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxcoef[n][s]*normal[n];
      for (m = 0; m < grid.dim; m++)
         bxcoef[m][s] = thesign*(ethere-ehere)/ethere*value*dotDndot[m]; //bxcoef[6][s] is coef for Du[s]
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxcoef[n][s]*Dndott[n];
      bxcoef[grid.dim][s] = thesign*(ethere-ehere)/ethere*value;
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxcoef[n][s]*Dndots[n];
      bxcoef[grid.dim+1][s] = thesign*(ethere-ehere)/ethere*value;
      bxcoef[grid.dim+2][s] = 0.0;
//      for (m = 0; m < 2*grid.dim; m++)
//         tempval[m] += bxcoef[m][s]*ux[s];
//         tempval[m] += bxcoef[m][s]*getDu(index,s,0,0,0.0,thesign,grid);

      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s][s]*normal[n];
      for (m = 0; m < grid.dim; m++)
         bxxcoef[m][s] = thesign*(ethere-ehere)/ethere*value*dotDndot[m]; //bxxcoef[6][s] is coef for DDu[s]
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s][s]*Dndott[n];
      bxxcoef[grid.dim][s] = thesign*(ethere-ehere)/ethere*
                             (value+normal[s]*tangent1[s]);
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s][s]*Dndots[n];
      bxxcoef[grid.dim+1][s] = thesign*(ethere-ehere)/ethere*
                               (value+normal[s]*tangent2[s]);
      bxxcoef[grid.dim+2][s] = 0.0;
//      for (m = 0; m < 2*grid.dim; m++)
//         tempval[m] += bxxcoef[m][s]*getD2u(index,s,s,0,0,0.0,thesign,grid);
   }

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 5)
   {
      for (m = 0; m < 2*grid.dim; m++)
         setvalarray(bcoef[m],sindex,0.0);
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            value = 0.0;
            for (s = 0; s < grid.dim; s++)
               value += D1uxxcoef[s][m][n]*normal[s];
            for (s = 0; s < grid.dim; s++)
               setvalarray(bcoef[s],sindex,
                           evalarray(bcoef[s],sindex)+
                           thesign*(ethere-ehere)/ethere*value*
                           evalarray(D2[m][n],sindex)*dotDndot[s]);
            value = 0.0;
            for (s = 0; s < grid.dim; s++)
               value += D1uxxcoef[s][m][n]*Dndott[s];
            setvalarray(bcoef[grid.dim],sindex,
                        evalarray(bcoef[grid.dim],sindex)+
                        thesign*(ethere-ehere)/ethere*
                        (value+normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                        evalarray(D2[m][n],sindex));
            value = 0.0;
            for (s = 0; s < grid.dim; s++)
               value += D1uxxcoef[s][m][n]*Dndots[s];
            setvalarray(bcoef[grid.dim+1],sindex,
                        evalarray(bcoef[grid.dim+1],sindex)+
                        thesign*(ethere-ehere)/ethere*
                        (value+normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                        evalarray(D2[m][n],sindex));
         }
//      for (s = 0; s < grid.dim; s++)
//         tindex[s] = index[s]-2+sindex[s];
//      for (s = 0; s < 2*grid.dim; s++)
//         tempval[s] += evalarray(bcoef[s],sindex)*getu(tindex,0,0,0.0,thesign,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 2;
//   cout << "b val ";
//   for (s = 0; s < 2*grid.dim; s++)
//      cout << tempval[s] << " ";
//   cout << endl;
// form part of d0 and dcoef's rth entry 
   gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);

   forwardbacksub0(temp,b0,LU,PLR,PLC,2*grid.dim-1);
   for (m = 0; m < grid.dim; m++)
      for (n = m; n < grid.dim; n++)
         u0[m][n] = temp[two2one[m][n]];
//   value = u0[0][0];
   for (s = 0; s < grid.dim; s++)
   {
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = bxcoef[m][s];
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      for (m = 0; m < grid.dim; m++)
         for (n = m; n < grid.dim; n++)
            uxcoef[m][n][s] = temp[two2one[m][n]];
//      value += uxcoef[0][0][s]*ux[s];
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = bxxcoef[m][s];
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      for (m = 0; m < grid.dim; m++)
         for (n = m; n < grid.dim; n++)
            uxxcoef[m][n][s] = temp[two2one[m][n]];
//      value += uxxcoef[0][0][s]*getD2u(index,s,s,0,0,0.0,thesign,grid);
   }

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 5)
   {
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = evalarray(bcoef[m],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      for (m = 0; m < grid.dim; m++)
         for (n = m; n < grid.dim; n++)
            setvalarray(ucoef[m][n],sindex,temp[two2one[m][n]]);
//      for (m = 0; m < grid.dim; m++)
//         tindex[m] = index[m]-2+sindex[m];
//      value += evalarray(ucoef[0][0],sindex)*getu(tindex,0,0,0.0,thesign,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 2;
//   cout << "value = " << value << endl;

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(bcoef[r],4,4,4);
}

void getjumpuxx(double &u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                int *index, int rstar, int sk, double alpha, int thesign, 
                double *normal, int mid, double ****D1ucoef, double **D1uxcoef, 
                double **D1uxxcoef, double*** D2[][3], double *D2uxcoef[][3], 
                double *D2uxxcoef[][3], double ***S, PBData &pb, GridData &grid)
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

   for (n = 0; n < grid.dim; n++)
      b0[n] = (sigma/ethere-Dtaudot[0])*dotDndot[n]+dotD2taudot[n];
   b0[grid.dim] = getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                  Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2];
   b0[grid.dim+1] = getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                    Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1];
   b0[grid.dim+2] = -jumpfe;

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
                         evalarray(D2[m][n],sindex));
            setvalarray(bcoef[grid.dim+1],sindex,
                        evalarray(bcoef[grid.dim+1],sindex)+
                        thesign*(ethere-ehere)/ethere*
                        (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                         evalarray(D2[m][n],sindex));
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

   gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);

   forwardbacksub0(temp,b0,LU,PLR,PLC,2*grid.dim-1);
   u0 = temp[two2one[rstar][rstar]];

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = evalarray(bcoef[m],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      setvalarray(ucoef,sindex,temp[two2one[rstar][rstar]]);

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
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = bxxcoef[m][s];
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      uxxcoef[s] = temp[two2one[rstar][rstar]];
   }

//   cout << "   in jumpuxx: " 
//        << evalcoef(b0,bcoef,bxcoef,bxxcoef,index,0,0,0.0,mid,thesign,grid) << endl;

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(bcoef[r],N,N,N);
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
// get Du in dim each dimension as linear combination of u and uxx, using D1 D2
// D1 is first order approxi of Du using u, D2 is approximation of cross-second-deriviative using u
// to calculate Du in dim r, uxxcoef[r][3] is coeff for [uxx, uyy, uzz], ucoef[r][3x3x3] is coeff for u
void getDu(double ****ucoef, double **uxxcoef, int *index, int rstar,
           int sstar, double alpha, int thesign, int *sk, double ****D1, 
           double*** D2[][3], GridData grid)
{
   int r, s, i;
   int sindex[grid.dim], tindex[grid.dim];

   for (r = 0; r < grid.dim; r++)
      for (s = 0; s < grid.dim; s++)
         uxxcoef[r][s] = 0.0;

   for (r = 0; r < grid.dim; r++)
   {
      uxxcoef[r][r] = 0.5*sk[r]*grid.dx[r];
      if (r == rstar)
         uxxcoef[r][r] += sk[r]*alpha*grid.dx[r];
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 0;
      while (sindex[0] < 3)
      {
         setvalarray(ucoef[r],sindex,evalarray(D1[r],sindex));
         if (r != rstar)
            setvalarray(ucoef[r],sindex,
                        evalarray(ucoef[r],sindex)+
                        sk[rstar]*alpha*grid.dx[rstar]*
                        evalarray(D2[min(r,rstar)][max(r,rstar)],sindex));

         (sindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
         {
            sindex[i] = 0;
            (sindex[i-1])++;
         }
      }
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 1;
   }

//   cout << "   check Du " << getDu(index,0,rstar,sstar,alpha,thesign,grid) << " " 
//        << getDu(index,1,rstar,sstar,alpha,thesign,grid) << " " 
//        << getDu(index,2,rstar,sstar,alpha,thesign,grid) << endl;
/*
   double D1val[grid.dim], D1more[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      D1val[r] = 0.0;
      D1more[r] = 0.0;
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 0;
      while (sindex[0] < 3)
      {
         for (s = 0; s < grid.dim; s++)
            tindex[s] = index[s]-1+sindex[s];
         D1val[r] += evalarray(D1[r],sindex)*getu(tindex,0,0,0.0,thesign,grid);
         D1more[r] += evalarray(ucoef[r],sindex)*getu(tindex,0,0,0.0,thesign,grid);
         (sindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
         {
            sindex[i] = 0;
            (sindex[i-1])++;
         }
      }
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 1;
      for (s = 0; s < grid.dim; s++)
         D1more[r] += uxxcoef[r][s]*getD2u(index,r,s,0,0,0.0,thesign,grid);
   }
*/
//   cout << "   check Du " << D1val[0] << " " << D1val[1] << " " << D1val[2] << endl;
//   cout << "   check Du " << D1more[0] << " " << D1more[1] << " " << D1more[2] << endl;
}
// In cim4, get Du-[r] as uxcoeff[r][3] ux uy uz and uxxcoeff[r][m][n] uxx uyy uzz
void getDu(double **uxcoef, double ***uxxcoef, int *index, int rstar, int sstar, 
           double alpha, int thesign, GridData grid)
{
   int r, s, m, i;
   int sindex[grid.dim], tindex[grid.dim];

   for (r = 0; r < grid.dim; r++)
      for (s = 0; s < grid.dim; s++)
      {
         uxcoef[r][s] = 0.0;
         for (m = 0; m < grid.dim; m++)
            uxxcoef[r][s][m] = 0.0;
      }

//   uxxcoef[rstar][rstar][rstar] = sstar*alpha*grid.dx[rstar];
   for (s = 0; s < grid.dim; s++)
   {
      uxcoef[s][s] = 1.0;
      uxxcoef[s][min(rstar,s)][max(rstar,s)] = sstar*alpha*grid.dx[rstar];
   }
}

void getDu(double *Du, int *index, int rstar, int sstar, double alpha, int thesign, 
           int *sk, double ****D1, double*** D2[][3], GridData grid)
{
   int r, s, i;
   int sindex[grid.dim], tindex[grid.dim];

   for (r = 0; r < grid.dim; r++)
   {
      Du[r] = 0.0;
      Du[r] += 0.5*sk[r]*grid.dx[r]*getD2u(index,r,r,0,0,0.0,thesign,grid);
      if (r == rstar)
         Du[r] += sk[r]*alpha*grid.dx[r]*getD2u(index,r,r,0,0,0.0,thesign,grid);
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 0;
      while (sindex[0] < 3)
      {
         for (s = 0; s < grid.dim; s++)
            tindex[s] = index[s]-1+sindex[s];
         Du[r] += evalarray(D1[r],sindex)*getu(tindex,0,0,0.0,thesign,grid);
         if (r != rstar)
            Du[r] += sk[rstar]*alpha*grid.dx[rstar]*
                     evalarray(D2[min(r,rstar)][max(r,rstar)],sindex)*
                     getu(tindex,0,0,0.0,thesign,grid);

         (sindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
         {
            sindex[i] = 0;
            (sindex[i-1])++;
         }
      }
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 1;
   }
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
      for (m = 2; m >= 1 && !foundone; m--)
         for (n = 2; n >= 1 && !foundone; n--)
            for (tempsk2[0] = -1,tempsk2[1] = tempsk2[0]+m; tempsk2[1] <= 1; 
                 (tempsk2[0])++,tempsk2[1] = tempsk2[0]+m)
               for (tempsk2[2] = -1,tempsk2[3] = tempsk2[2]+n; tempsk2[3] <= 1; 
                    (tempsk2[2])++,tempsk2[3] = tempsk2[2]+n)
               {
// check the four points of tempsk2
                  bad = 0;
                  for (r = 0; r < 2 && !bad; r++)
                  {
                     rindex[i] = index[i]+tempsk2[r];
                     for (s = 2; s < 4 && !bad; s++)
                     {
                        rindex[j] = index[j]+tempsk2[s];
                        if (rindex[i] < 0 || rindex[i] > grid.nx[i] || 
                            rindex[j] < 0 || rindex[j] > grid.nx[j] || 
                            (evalarray(S,index) < 0.0)+
                            (evalarray(S,rindex) < 0.0) == 1)
                           bad = 1;
                        else// if for current combimation of tempsk2, rindex is not bad(change sign), record status of rindex, sort from large to small
                        {// that is, tempsk2stat[0] > ... tempsk2stat[3]
// getstatus uses yessk2
                           tempsk2stat[2*r+(s-2)] = getstatus5(S,rindex,grid);
                           for (t = 2*r+(s-2); t > 0 && 
                                               tempsk2stat[t] > tempsk2stat[t-1]; t--)
                           {//goes down the list, swap if tempsk2stat[t]>tempsk2stat[t-1]
                              temp = tempsk2stat[t-1];
                              tempsk2stat[t-1] = tempsk2stat[t];
                              tempsk2stat[t] = temp;
                           }
                        }
                     }
                     rindex[j] = index[j];
                  }
                  rindex[i] = index[i];
               
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

#if 0
char yessk2(int *sk2, int i, int j, int *index, double ***S, GridData &grid)
{
   int r, s, t, m, n, temp, rindex[grid.dim], tempsk2[4], sk2stat[4], tempsk2stat[4];
   char bad;

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];

   bad = 1;
   for (m = 2; m >= 1 && bad; m--)
      for (tempsk2[0] = -1,tempsk2[1] = tempsk2[0]+m; tempsk2[1] <= 1 && bad; 
           (tempsk2[0])++,tempsk2[1] = tempsk2[0]+m)
         for (n = 2; n >= 1 && bad; n--)
            for (tempsk2[2] = -1,tempsk2[3] = tempsk2[2]+n; tempsk2[3] <= 1 && bad; 
                 (tempsk2[2])++,tempsk2[3] = tempsk2[2]+n)
            {
// check the four points of tempsk2
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
               
               if (!bad)
                  for (r = 0; r < 4; r++)
                     sk2[r] = tempsk2[r];
            }

   return (!bad);
}
#endif
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

int getstatus4(double ***S, int *index, GridData grid)
{
   int r, s, tempsk2[4];
   int thestatus, except;
   
   thestatus = getstatus3(S,index,grid);
   if (thestatus == 2 || thestatus == 3)
   {
      except = 0;
      for (r = 0; r < grid.dim && !except; r++)
         for (s = r+1; s < grid.dim && !except; s++)
            if (!yessk2(tempsk2,r,s,index,S,grid))
               except = 1; //if not sk2, then is exacpt, status = 3
      if (!except)
         thestatus = 2;
      else
         thestatus = 3;
   }

   return thestatus;
}
//cim3 method to get interface gradient
double getinterfacegrad3(double *grad, double ***u, double ***S, int *index, int rstar, 
                         int sstar, PBData &pb, GridData grid)
{
   int r, s, i, j, k, tmps;
   int rindex[grid.dim], tindex[grid.dim];
   double **LU, **M, value, Sval; 
   int PLR[grid.dim], PLC[grid.dim];
   double alpha, beta, tangent[grid.dim], normal[grid.dim], temp[grid.dim], 
          a[4], m[grid.dim];
   double bk, rthere, rhere, ethere, ehere, ehat, ebar, C;
   int thestatus, gamma[grid.dim][2], sk[grid.dim];
   double ouralpha;

   for (r = 0; r < grid.dim; r++)
      grad[r] = 0.0;

   Sval = evalarray(S,index);
   for (r = 0; r < grid.dim; r++)
      tindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
   {
      for (s = -1; s <= 1; s += 2)
      {
         tindex[r] = min(max(index[r]+s,0),grid.nx[r]);
         if (Sval*evalarray(S,tindex) < 0.0)
            gamma[r][(s+1)/2] = 1;
         else
            gamma[r][(s+1)/2] = 0;
      }
      sk[r] = gamma[r][1]-gamma[r][0];
      tindex[r] = index[r];
   }

   if (gamma[rstar][(sstar+1)/2] == 0)
      return 0;

   LU = matrix(grid.dim-1,grid.dim-1);
   M = matrix(grid.dim-1,grid.dim-1);

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

   thestatus = getstatus3(S,index,grid);
   if (thestatus == 2)
   {
      for (r = 0; r < grid.dim; r++)
         temp[r] = 0.0;
      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];
      for (r = 0; r < grid.dim; r++)
      {
         rindex[r] = index[r]+sk[r];
         getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
         if (r == rstar)
            ouralpha = alpha;
         rindex[r] = index[r];
         if (abs(sk[r]) == 1)
         {
            beta = 1.0-alpha;
            ehat = (beta+beta*beta)*(0.5+alpha)*ehere+
                   (alpha+alpha*alpha)*(0.5+beta)*ethere;
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
               temp[r] += a[s]*evalarray(u,rindex);
            }
            rindex[r] = index[r]-sk[r];
            temp[r] += -sk[r]*C*tangent[r]*sk[r]*evalarray(u,rindex);
            if (r == rstar)
               grad[r] += -sk[r]*evalarray(u,rindex);
            rindex[r] = index[r];
            temp[r] += sk[r]*C*tangent[r]*sk[r]*evalarray(u,rindex);
            if (r == rstar)
               grad[r] += sk[r]*evalarray(u,rindex);
            M[r][r] = 1.0-abs(sk[r])*(0.5+alpha)*bk*tangent[r]*tangent[r];
            if (r == rstar)
               m[r] = sk[r]*(0.5+alpha);
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
                     m[j] = 0.5*tmps;

                  if (abs(tmps) == 1)
                  {
                     rindex[j] = index[j]-tmps;
                     temp[r] += -sk[r]*C*tangent[j]*(1.0+alpha)*tmps*
                                 evalarray(u,rindex);
                     if (r == rstar)
                        grad[j] += -tmps*(1.0+alpha)*evalarray(u,rindex);
                     rindex[r] = index[r]-sk[r];
                     temp[r] += sk[r]*C*tangent[j]*alpha*tmps*evalarray(u,rindex);
                     if (r == rstar)
                        grad[j] += tmps*alpha*evalarray(u,rindex);
                     rindex[j] = index[j];
                     temp[r] += -sk[r]*C*tangent[j]*alpha*tmps*evalarray(u,rindex);
                     if (r == rstar)
                        grad[j] += -tmps*alpha*evalarray(u,rindex);
                     rindex[r] = index[r];
                     temp[r] += sk[r]*C*tangent[j]*(1.0+alpha)*tmps*
                                evalarray(u,rindex);
                     if (r == rstar)
                        grad[j] += tmps*(1.0+alpha)*evalarray(u,rindex);
                  }
                  else
                  {
                     for (s = -1; s <= 1; s += 2)
                     {
                        rindex[j] = index[j]+s;
                        temp[r] += sk[r]*C*tangent[j]*(1.0+alpha)*0.5*s*
                                   evalarray(u,rindex);
                        if (r == rstar)
                           grad[j] += s*(1.0+alpha)*0.5*evalarray(u,rindex);
                     }
                     rindex[r] = index[r]-sk[r];
                     for (s = -1; s <= 1; s += 2)
                     {
                        rindex[j] = index[j]+s;
                        temp[r] += -sk[r]*C*tangent[j]*alpha*0.5*s*evalarray(u,rindex);
                        if (r == rstar)
                           grad[j] += -s*alpha*0.5*evalarray(u,rindex);
                     }
                     rindex[r] = index[r];
                     rindex[j] = index[j];
                  }
               }
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
               temp[r] += evalarray(u,rindex);
            }
            rindex[r] = index[r];
            temp[r] += -2.0*evalarray(u,rindex);
         }
      }

      gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
   
      for (r = 0; r < grid.dim; r++)
         grad[r] += m[r]*temp[r];
      for (r = 0; r < grid.dim; r++)
         grad[r] /= grid.dx[rstar];

      rindex[rstar] = index[rstar]-sstar;
      value = evalarray(u,index)+
              ((evalarray(u,index)-evalarray(u,rindex))/grid.dx[rstar]+
               0.5*grid.dx[rstar]*temp[rstar])*ouralpha*grid.dx[rstar]+
              0.5*temp[rstar]*ouralpha*ouralpha*grid.dx[rstar]*grid.dx[rstar];
   }
   else if (thestatus == 3)
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];
      s = sstar;
      for (r = 0; r < grid.dim; r++)
      {
         rindex[r] = index[r]+s;
         getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
         if (r == rstar)
            ouralpha = alpha;
         rindex[r] = index[r];
         beta = 1.0-alpha;
         ebar = alpha*ethere+beta*ehere;
         for (j = 0; j < grid.dim; j++)
            M[r][j] = -gamma[r][(s+1)/2]*gamma[j][(s+1)/2]*beta*(ehere-ethere)/ebar*
                       tangent[r]*tangent[j];
         M[r][r] += 1.0;
            
         for (j = 0; j < grid.dim; j++)
         {
            rindex[j] = index[j]+s;
            grad[r] += (gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                        (1-gamma[j][(s+1)/2])*s*tangent[j])/(grid.dx[r]*ebar)*
                       evalarray(u,rindex);
            rindex[j] = index[j];
         }
         rindex[r] = index[r]+s;
         grad[r] += (s*ethere)/(grid.dx[r]*ebar)*evalarray(u,rindex);
         rindex[r] = index[r];
         grad[r] += -(s*ethere)/(grid.dx[r]*ebar)*evalarray(u,rindex);
         for (j = 0; j < grid.dim; j++)
            grad[r] += -(gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                         (1-gamma[j][(s+1)/2])*s*tangent[j])/(grid.dx[r]*ebar)*
                        evalarray(u,rindex);
      }

      gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);
      forwardbacksub0(grad,grad,LU,PLR,PLC,grid.dim-1);

      value = evalarray(u,index)+grad[rstar]*sstar*ouralpha*grid.dx[rstar];
   }
   else
      cout << "Error" << endl;

   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);

   return value;
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

void linearsystembasic(SparseElt2**** &A, double ***b, double ***S, GridData &grid)
{
// can use cim1, cim2, cim3, cim4, cim5 and saves on speed in calculating Du
   int i, r, s, j, k, thestatus;
   int tindex[grid.dim], rindex[grid.dim], numcount = 3, count[numcount];
   double x[grid.dim];
   SparseElt2 *current;

   cout << "in linearsystembasic" << endl;

   clearsparse(A,grid);

   for (i = 0; i < numcount; i++)
      count[i] = 0;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(b,tindex,1.0);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      thestatus = 1;
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim && thestatus == 1; r++)
      {
         for (s = -1; s <= 1 && thestatus == 1; s += 2)
         {
            rindex[r] = tindex[r]+s;
            if (rindex[r] < 0 || rindex[r] > grid.nx[r])
               thestatus = 0;//if neighbor out of bound, tindex has status 0
         }
         rindex[r] = tindex[r];
      }
      if (thestatus == 1)//inbox
      {
         double tempval = 0.0;
         for (r = 0; r < grid.dim; r++)
            rindex[r] = tindex[r];
         for (r = 0; r < grid.dim && thestatus == 1; r++)
         {
            tempval += 2.0/grid.dx[r]*grid.dx[r];
            for (s = -1; s <= 1 && thestatus == 1; s += 2)
            {
               rindex[r] = tindex[r]+s;
               sparse2(tindex,rindex,A,-1.0/(grid.dx[r]*grid.dx[r]),grid);
            }
            rindex[r] = tindex[r];
         }
         sparse2(tindex,tindex,A,tempval,grid);
         (count[1])++;
      }
      else if (thestatus == 0)//boundary
      {
         sparse2(tindex,tindex,A,1.0,grid);
         for (r = 0; r < grid.dim; r++)
            x[r] = grid.a[r]+tindex[r]*grid.dx[r];
         setvalarray(b,tindex,1.0);
         (count[0])++;
      }
      else
      {
         cout << "Error in status" << endl;
         exit(1);
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "found " << count[0] << " boundary pts and " << count[1] << " interior pts."
        << endl;
}
// used to solve pb
void linearsystem3(SparseElt2**** &A, double ***b, double ***S, PBData &pb, 
                   GridData &grid)
{
// more flexible hybrid CIM2 designation
   int gamma[grid.dim][2];
   int i, r, s, j, k, thestatus;
   int tindex[grid.dim], rindex[grid.dim], count[3];
   double Sval1, Sval2;
   double x[grid.dim];
   char except;
   double temp;

   clearsparse(A,grid);

   for (i = 0; i < grid.dim; i++)
      count[i] = 0;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim && tindex[r] > 0 && tindex[r] < grid.nx[r]; r++);
      if (evalarray(S,tindex) < 0.0)
         setvalarray(b,tindex,0.0);
      else
      {
         if (r == grid.dim)
         {
            sub2coord(x,tindex,grid);
            temp = evalarray(pb.psi,tindex)+getpsivac(x,pb);
            setvalarray(b,tindex,-Bprime(temp,pb)+
                                  B2prime(temp,pb)*evalarray(pb.psi,tindex));
            sparse2(tindex,tindex,A,B2prime(temp,pb),grid);
         }
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

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
      thestatus = getstatus3(S,tindex,grid);

      if (thestatus == 1)
      {
         interiorpt2(A,b,tindex,S,pb,grid);
         (count[0])++;
      }
      else if (thestatus == 2)
      {
         cim2again3(A,b,tindex,gamma,S,pb,grid);
         (count[2])++;
      }
      else if (thestatus == 3)
      {
         cim1again2(A,b,tindex,gamma,S,pb,grid);
         (count[1])++;
      }
      else
      {
         sparse2(tindex,tindex,A,1.0,grid);
         for (r = 0; r < grid.dim; r++)
            x[r] = grid.a[r]+tindex[r]*grid.dx[r];
         setvalarray(b,tindex,DBC(x,grid.dim,0.0));
         (count[0])++;
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "found " << count[0] << " interior pts, " << count[1] << " cim1 pts, "
        << count[2] << " cim2 pts." << endl;
}
// No Poisson Boltzman, no Dusmall
void linearsystem4(SparseElt2**** &A, double ***b, double ***S, PBData &pb, 
                   GridData &grid)
{
// can use cim1, cim2, cim3, cim4, cim5
   int gamma[grid.dim][2];
   int i, r, s, j, k, thestatus;
   int tindex[grid.dim], rindex[grid.dim], numcount = 7, count[numcount];
   double Sval1, Sval2;
   double x[grid.dim];
   char except;
   double temp;
   int tempsk2[4];

   cout << "in linearsystem4" << endl;

   clearsparse(A,grid);

   for (i = 0; i < numcount; i++)
      count[i] = 0;

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
      if (globcim == 2)
         thestatus = getstatus3(S,tindex,grid);
      else if (globcim == 4)
         thestatus = checkcim3(S,tindex,grid);
      else if (globcim == 5 || globcim == 345)
         thestatus = getstatus5(S,tindex,grid);

      if (thestatus == 1)
      {
         interiorptsmall(A,b,tindex,S,pb,grid);
         (count[1])++;
      }
      else if (thestatus == 2)
      {
         if (globcim == 2)
            cim2again3(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 4 || globcim == 5 || globcim == 345)
         {
            if (!globallcim345)
               cim3again(A,b,tindex,gamma,S,pb,grid);
            else
               cim345(A,b,tindex,gamma,S,pb,grid);
         }
         (count[2])++;
      }
      else if (thestatus == 3)
      {
         if (globcim == 2)
            cim1again2(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 4)
            cim4(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 5)
            cim5(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 345)
            cim345(A,b,tindex,gamma,S,pb,grid);
         (count[3])++;
      }
      else if (thestatus == 4)
      {
         if (globcim == 4)
            cim1again2(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 5)
         {
/*
            for (j = -1; j <= 1; j++)
            {
               for (s = 1; s >= -1; s--)
               {
                  for (r = -1; r <= 1; r++)
                     cout << S[tindex[0]+r][tindex[1]+s][tindex[2]+j] << " ";
                  cout << endl;
               }
               cout << endl;
            }
            exit(1);
*/
            cim4(A,b,tindex,gamma,S,pb,grid);
         }
         else if (globcim == 345)
            cim345(A,b,tindex,gamma,S,pb,grid);
         (count[4])++;
      }
      else if (thestatus == 5)
      {
         cim1again2(A,b,tindex,gamma,S,pb,grid);
         (count[5])++;
      }
      else if (thestatus == 0)//boundary point
      {
         sparse2(tindex,tindex,A,1.0,grid);
         for (r = 0; r < grid.dim; r++)
            x[r] = grid.a[r]+tindex[r]*grid.dx[r];
         setvalarray(b,tindex,DBC(x,grid.dim,0.0));
         (count[0])++;
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
//          thecount2 += sizeof(SparseElt2)+sizeof(int)*3;
      }
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
//   cout << "number of matrix entries = " << thecount << " for " << thecount*36 << endl;
   cout << "found " << count[0] << " boundary pts, " << count[1] << " interior pts, ";
   if (globcim == 2)
      cout << count[2] << " cim2 pts, " << count[3] << " cim1 pts." << endl;
   else if (globcim == 4)
      cout << count[2] << " cim3 pts, " << count[3] << " cim4 pts," << endl
           << "   and " << count[4] << " cim1 pts." << endl;
   else if (globcim == 5 || globcim == 345)
      cout << count[2] << " cim3 pts, " << count[3] << " cim5 pts," << endl
           << "   and " << count[4] << " cim4 pts, " << count[5] << " cim1 pts." 
           << endl;
/*
   ofstream chi2d("chi2d.dat", ios::out);
   outputsparsesmall(chi2d,A,S,pb,grid);
   ofstream gamma2d("gamma2d.dat", ios::out);
   gamma2d.precision(16);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      gamma2d << scientific << evalarray(b,tindex) << endl;
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
*/
}
// same as linearsystem4, use StorageStruct
void linearsystem5(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, 
                   int &smallsize, double ***S, PBData &pb, GridData &grid)
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

   cout << "in linearsystem5" << endl;

   clearsparse(A,grid);

   int buildsize = 0;
   smallsize = getstoragesize(S,grid);
   Dusmall = new StorageStruct[smallsize];
//   for (i = 0; i < smallsize; i++)
//   {
//      Dusmall[i].N = 4;
//      Dusmall[i].info = new int[Dusmall[i].N];
//      Dusmall[i].head = NULL;
//   }
   cout << "smallsize = " << smallsize << endl;

   for (i = 0; i < numcount; i++)
      count[i] = 0;

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
            rindex[r] = min(max(tindex[r]+s,0),grid.nx[r]);//if rindex (nbr of tindex) goes out of boundary, set as tindex
            if (evalarray(S,tindex)*evalarray(S,rindex) < 0.0)
               gamma[r][(s+1)/2] = 1;
            else
               gamma[r][(s+1)/2] = 0;
         }
         rindex[r] = tindex[r];
      }
      if (globcim == 2)
         thestatus = getstatus3(S,tindex,grid);
      else if (globcim == 4)
         thestatus = checkcim3(S,tindex,grid);
      else if (globcim == 5 || globcim == 345|| globcim == 6)
         thestatus = getstatus5(S,tindex,grid);

      if (thestatus == 1)
      {
         interiorptsmall(A,b,tindex,S,pb,grid);
         (count[1])++;
      }
      else if (thestatus == 2)
      {
         if (globcim == 2)
            cim2(A,b,Dusmall,buildsize,tindex,gamma,S,pb,grid);
//            cim2again3(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 4 || globcim == 5 || globcim == 345)
         {
            if (!globallcim345)
               cim3again(A,b,tindex,gamma,S,pb,grid);
            else
               cim345(A,b,Dusmall,buildsize,tindex,gamma,S,pb,grid);
         }
         else if (globcim == 6)
            cim345cond(A,b,Dusmall,buildsize,tindex, nullptr, gamma,S,pb,grid);
         else if (globcim == 1)
            cim1(A,b,Dusmall,buildsize,tindex,gamma,S,pb,grid);
         (count[2])++;
      }
      else if (thestatus == 3)
      {
         if (globcim == 2 || globcim == 1)
            cim1(A,b,Dusmall,buildsize,tindex,gamma,S,pb,grid);
//            cim1again2(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 4)
            cim4(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 5)
            cim5(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 345)
            cim345(A,b,Dusmall,buildsize,tindex,gamma,S,pb,grid);
          else if (globcim == 6)
            cim345cond(A,b,Dusmall,buildsize,tindex, nullptr, gamma,S,pb,grid);
         (count[3])++;
      }
      else if (thestatus == 4)
      {
         if (globcim == 4)
            cim1(A,b,Dusmall,buildsize,tindex,gamma,S,pb,grid);
//            cim1again2(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 5)
            cim4(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 345)
            cim345(A,b,Dusmall,buildsize,tindex,gamma,S,pb,grid);
          else if (globcim == 6)
            cim345cond(A,b,Dusmall,buildsize,tindex, nullptr, gamma,S,pb,grid);
         (count[4])++;
      }
      else if (thestatus == 5)
      {
         cim1(A,b,Dusmall,buildsize,tindex,gamma,S,pb,grid);
//         cim1again2(A,b,tindex,gamma,S,pb,grid);
         (count[5])++;
      }
      else if (thestatus == 0) //boundary, set A diag = 1, set b = boundary
      {
         sparse2(tindex,tindex,A,1.0,grid);
         for (r = 0; r < grid.dim; r++)
            x[r] = grid.a[r]+tindex[r]*grid.dx[r];
         setvalarray(b,tindex,DBC(x,grid.dim,0.0));
         (count[0])++;
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
   if (globcim == 2)
      cout << count[2] << " cim2 pts, " << count[3] << " cim1 pts." << endl;
   else if (globcim == 4)
      cout << count[2] << " cim3 pts, " << count[3] << " cim4 pts," << endl
           << "   and " << count[4] << " cim1 pts." << endl;
   else if (globcim == 5 || globcim == 345|| globcim == 6)
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

// USING EXACT
//   double ***res = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
//   cout << "exact residual = " << getexactresidual(res,A,b,grid,S,pb) << endl;
//   output(mu2d,res,grid.nx);
//   free_matrix(res,grid.nx[0],grid.nx[1],grid.nx[2]);
//   exit(1);

   // output linear system
   // ofstream chi2d("chi2d.dat", ios::out);
   // outputsparsesmall(chi2d,A,S,pb,grid);
   // ofstream gamma2d("gamma2d.dat", ios::out);
   // gamma2d.precision(16);
   // for (i = 0; i < grid.dim; i++)
   //    tindex[i] = 0;
   // while (tindex[0] <= grid.nx[0])
   // {
   //    gamma2d << scientific << evalarray(b,tindex) << endl;
   //    (tindex[grid.dim-1])++;
   //    for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
   //    {
   //       tindex[i] = 0;
   //       (tindex[i-1])++;
   //    }
   // }
}
// with a, use dusmall
void linearsystem6(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, 
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

   cout << "in linearsystem6" << endl;

   clearsparse(A,grid);

   int buildsize = 0;
   smallsize = getstoragesize(S,grid);
   Dusmall = new StorageStruct[smallsize];
//   for (i = 0; i < smallsize; i++)
//   {
//      Dusmall[i].N = 4;
//      Dusmall[i].info = new int[Dusmall[i].N];
//      Dusmall[i].head = NULL;
//   }
   cout << "smallsize = " << smallsize << endl;

   for (i = 0; i < numcount; i++)
      count[i] = 0;

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
      if (globcim == 2)
         thestatus = getstatus3(S,tindex,grid);
      else if (globcim == 4)
         thestatus = checkcim3(S,tindex,grid);
      else if (globcim == 5 || globcim == 345 || globcim == 6)
         thestatus = getstatus5(S,tindex,grid);

      if (thestatus == 1)
      {
         interiorptsmall(A,b,tindex,S,pb,grid);
         (count[1])++;
      }
      else if (thestatus == 2)
      {
         if (globcim == 2)
            cim2(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
         else if (globcim == 4 || globcim == 5 || globcim == 345)
         {
            if (!globallcim345)
               cim3again(A,b,tindex,gamma,S,pb,grid);
            else
               cim345(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
         }
         else if (globcim == 6)
            cim345cond(A,b,Dusmall,buildsize,tindex, a, gamma,S,pb,grid);
         else if (globcim == 1)
            cim1(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
         (count[2])++;
      }
      else if (thestatus == 3)
      {
         if (globcim == 2 || globcim == 1)
            cim1(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
         else if (globcim == 4)
            cim4(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 5)
            cim5(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 345)
            cim345(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
          else if (globcim == 6)
            cim345cond(A,b,Dusmall,buildsize,tindex, a, gamma,S,pb,grid);
         (count[3])++;
      }
      else if (thestatus == 4)
      {
         if (globcim == 4)
            cim1(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
//            cim1again2(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 5)
            cim4(A,b,tindex,gamma,S,pb,grid);
         else if (globcim == 345)
            cim345(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
          else if (globcim == 6)
            cim345cond(A,b,Dusmall,buildsize,tindex, a, gamma,S,pb,grid);
         (count[4])++;
      }
      else if (thestatus == 5)
      {
         cim1(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
//         cim1again2(A,b,tindex,gamma,S,pb,grid);
         (count[5])++;
      }
      else if (thestatus == 0)
      {
         sparse2(tindex,tindex,A,1.0,grid);
         for (r = 0; r < grid.dim; r++)
            x[r] = grid.a[r]+tindex[r]*grid.dx[r];
         setvalarray(b,tindex,DBC(x,grid.dim,0.0));
         (count[0])++;
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
   if (globcim == 2)
      cout << count[2] << " cim2 pts, " << count[3] << " cim1 pts." << endl;
   else if (globcim == 4)
      cout << count[2] << " cim3 pts, " << count[3] << " cim4 pts," << endl
           << "   and " << count[4] << " cim1 pts." << endl;
   else if (globcim == 5 || globcim == 345 || globcim == 6)
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

// USING EXACT
//   double ***res = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
//   cout << "exact residual = " << getexactresidual(res,A,b,grid,S,pb) << endl;
//   output(mu2d,res,grid.nx);
//   free_matrix(res,grid.nx[0],grid.nx[1],grid.nx[2]);
//   exit(1);
}

void linearsystemZL(SparseElt2**** &A, double ***b, double ***S, PBData &pb, 
                    GridData &grid)
{
// can use cim1, cim2, cim3, cim4, cim5 and saves on speed in calculating Du
   int i, r, s, j, k, thestatus;
   int tindex[grid.dim], rindex[grid.dim], numcount = 3, count[numcount];
   double x[grid.dim];
   SparseElt2 *current;

   cout << "in linearsystemZL" << endl;

   clearsparse(A,grid);

//   int buildsize = 0;
//   smallsize = getstoragesize(S,grid);
//   Dusmall = new StorageStruct[smallsize];
//   cout << "smallsize = " << smallsize << endl;

   for (i = 0; i < numcount; i++)
      count[i] = 0;

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

   char ***tube = cmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(tube,tindex,1);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

//   cout << "testing ZL" << endl;
//   testZL(0,5,S,tube,pb,grid);

   double tempval, maxtempval = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      thestatus = 1;
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim && thestatus == 1; r++)
      {
         for (s = -1; s <= 1 && thestatus == 1; s += 2)
         {
            rindex[r] = tindex[r]+s;
            if (rindex[r] >= 0 && rindex[r] <= grid.nx[r])
               if (evalarray(S,tindex)*evalarray(S,rindex) < 0.0)
                  thestatus = 2;
               else;
            else
               thestatus = 0;
         }
         rindex[r] = tindex[r];
      }
      if (thestatus == 1)
      {
         interiorptsmall(A,b,tindex,S,pb,grid);
         (count[1])++;
      }
      else if (thestatus == 2)
      {
//         if (tindex[0] == 10 && tindex[1] == 20 && tindex[2] == 20)
//         {
//            globdebug = 1;
//            cout << "STARTING DEBUG" << endl;
//         }
//         globdebug = 1;
//         iim(A,b,tindex,0,5,S,tube,pb,grid);
         iimghost(A,b,tindex,0,5,S,tube,pb,grid);
//         globdebug = 0;
//         getchar();
         (count[2])++;
      }
      else if (thestatus == 0)
      {
         sparse2(tindex,tindex,A,1.0,grid);
         for (r = 0; r < grid.dim; r++)
            x[r] = grid.a[r]+tindex[r]*grid.dx[r];
         setvalarray(b,tindex,DBC(x,grid.dim,0.0));
         (count[0])++;
      }
      else
      {
         cout << "Error in status" << endl;
         exit(1);
      }

      if (thestatus == 2)
      {
         tempval = evalarray(b,tindex);
         for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
            if (evalarray(S,(*current).cindex) < 0.0)
               tempval -= (*current).val*getu((*current).cindex,0,0,0.0,-1,grid);
            else
               tempval -= (*current).val*getu((*current).cindex,0,0,0.0,1,grid);
         if (fabs(tempval) > maxtempval)
            maxtempval = fabs(tempval);
//         if (maxtempval > 80.0)
//         {
//            cout << tindex[0] << " " << tindex[1] << " " << tindex[2] << endl;
//            getchar();
//         }
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "RESIDUAL FROM EXACT = " << maxtempval << endl;

   cout << "found " << count[0] << " boundary pts, " << count[1] << " interior pts, ";
   cout << count[2] << " interface pts." << endl;

   free_matrix(tube,grid.nx[0],grid.nx[1],grid.nx[2]);
// Dusmall will hold more memory: smallsize; but actually use less: buildsize
//   smallsize = buildsize;
//   cout << "smallsize reduced to " << smallsize << endl;

//   thecount2 = 0;
//   for (i = 0; i < smallsize; i++)
//      for (SparseElt *current = Dusmall[i].head; current != NULL; 
//           current = (*current).next)
//         thecount2++;
//   cout << "Du storage structure has " << thecount2 << " elements." << endl;

// USING EXACT
//   double ***res = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
//   cout << "exact residual = " << getexactresidual(res,A,b,grid,S,pb) << endl;
//   output(mu2d,res,grid.nx);
//   free_matrix(res,grid.nx[0],grid.nx[1],grid.nx[2]);
//   exit(1);
}
// 0=boundary, 1=interior, 2= can cim2, 3 = need cim3, 4 = need cim4
int checkcim3(double ***S, int *index, GridData &grid)
{
   int r, sk, m, n, rindex[grid.dim], sk2[4];
   char yescim3 = 1, yesinterior = 1;

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
      {
         rindex[r] = index[r]+sk;
         if (rindex[r] < 0 || rindex[r] > grid.nx[r])
            return 0;
         rindex[r] = index[r];
      }

   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
      {
         rindex[r] = index[r]+sk;
         if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)
            yesinterior = 0;
         rindex[r] = index[r];
      }
   if (yesinterior)
      return 1;

   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
         if (!yessk2(sk2,m,n,index,S,grid)) // yessk2 is false
         {
            yescim3 = 0;
            for (r = 0; r < grid.dim; r++)
               for (sk = -1; sk <= 1; sk += 2)
               {
                  rindex[r] = index[r]+sk;
                  if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1 &&
                      !yessk2(sk2,m,n,rindex,S,grid))
                     return 4;//look at points across interface, but that point has yessk2=false
                  rindex[r] = index[r];
               }
         }
   if (yescim3)// yessk2 is true for every plane
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
         return 2;//one interface point, can use cim2,
   }
   return 3;// look at points across interface, those points all have yessk2=true, or count >= 2
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
//cim12.pdf p35, approximate mixed derivtive in (m,n) plane at index, 
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
   thesign = 2*(evalarray(S,tindex) >= 0.0)-1;//1=outside, -1 indside
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

// looks for central differencing possibility around node
   for (r = 0; r < ntheta; r += 2)
      if ((thesign < 0)+(signs[r] < 0) != 1 && 
          (thesign < 0)+(signs[(r+2)%ntheta] < 0) != 1)
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

// looks for forward differencing possibility around node
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
               int count = 0;
               for (r = 0; r < grid.dim; r++)
                  for (sk = -1; sk <= 1; sk += 2)
                  {
                     rindex[r] = index[r]+sk;
                     if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1 &&
                         !yessk2(sk2,m,n,rindex,S,grid)){
                          ++count;
                          // cimstatus[m][n] = 0;
                      }
                     rindex[r] = index[r];
                  }
                if (count==6){
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

char checkcimstatus(double ***S, GridData &grid)
{
   int i, temp1, temp2, tindex[grid.dim];

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp1 = getstatus4(S,tindex,grid);
      temp2 = checkcim3(S,tindex,grid);
      if (temp2 == 4 || ((temp1 <= 1 || temp2 <= 1) && temp1 != temp2) || 
          ((temp1 > 1 || temp2 > 1) && temp1 != temp2))
         return 0;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   return 1;
}


//get rhs of level set equation. without a term. output is rhs, u is surface,
void getRHS(double ***rhs, double ***u, PBData &pb, MarchStruct &march, TempStruct &tmp, 
            GridData &grid)
{
   double dfp[grid.dim], dfn[grid.dim], df[grid.dim];
   double grad, maxvn = 0.0, exactmaxvn;
   int deg = 3;
   double u1d[2*deg+1];
   int i, r, s, j, k;
   int tindex[grid.dim], rindex[grid.dim];
   double Du[grid.dim], D2u[grid.dim][grid.dim];
   double value, grad2, tol = grid.tol;
   // ofstream phi2d("phi2d.dat", ios::out);
   // ofstream beta2d("beta2d.dat", ios::out);
   clock_t cstart, cend;
   double*** sign_surf;// For ICIM
   for (r = 0; r < grid.dim; r++)
   {
      globerrvec[r] = 0.0;
      globerrvec2[r] = 0.0;
      globerrvec3[r] = 0.0;
      globerrvec4[r] = 0.0;
      for (s = 0; s < grid.dim; s++)
         globerrmx[r][s] = 0.0;
   }

   cout << "in getRHS" << endl;
// CHANGE HERE
// solves the Poisson-Boltzman equation and puts the answer in the 3D array pb.psi
//   linearsystem3(tmp.A,tmp.b,u,pb,grid);
   cstart = clock ();
//   if (checkcimstatus(u,grid))
//      cout << "checks out" << endl;
//   else
//   {
//      cout << "bad checking" << endl;
//      exit(1);
//   }
   StorageStruct *Dusmall;
   int smallsize = 0;
   if (!globsmall)
      if (globcim != 0)//cim
         linearsystem4(tmp.A,tmp.b,u,pb,grid);//interior point not central difference, put in matrix solver
          // linearsystem3(tmp.A,tmp.b,u,pb,grid);
      else//iim
      {
         // linearsystemZL(tmp.A,tmp.b,u,pb,grid);
        linearsystembasic(tmp.A,tmp.b,u,grid);
      }
   else if (globcim !=7) {
      linearsystem5(tmp.A,tmp.b,Dusmall,smallsize,u,pb,grid);
      march.Dusmall = Dusmall;
      march.smallsize = smallsize;
   }else if (globcim == 7){
      sign_surf = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
      SignSurf(sign_surf, u, grid);
      int maxIter = 10;
      int depth = 1;
      bool flipbydet = true;
      flip(sign_surf, u, pb, grid, maxIter, depth, flipbydet);
      linearsystem_icim(tmp.A,tmp.b,u, sign_surf, pb, grid);
   }else{
      cerr<<"No available option for getRHS" <<endl;
      exit(1);
   }
/*
   ofstream gamma2d("gamma2d.dat", ios::out);
   for (r = 0; r < smallsize; r++)
   {
      gamma2d << Dusmall[r].info[0] << " " << Dusmall[r].info[1] << " "
              << Dusmall[r].info[2] << " " << Dusmall[r].info[3] << endl;
   }
   gamma2d.close();
*/
   cend = clock ();
   cout << "linear system clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC 
        << endl;
   cout << "equation error = " << globerr << " at " << globerrvec3[0] << " " 
        << globerrvec3[1] << " " << globerrvec3[2] << endl;
//   BICGSTAB(pb.psi,tmp.A,tmp.b,grid,grid.nx[0]*grid.nx[1]*grid.nx[2],1.0e-14,tmp);
/*
   SparseElt2 ****B; 
   B = sparseelt2ptrmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   for (r = 0; r < (grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1); r++)
   {
      ind2sub(rindex,r,grid.nx,grid.dim);
      setvalarray(tmp.b,rindex,1.0);
      if (rindex[0] == 0 || rindex[1] == 0 || rindex[2] == 0 || 
          rindex[0] == grid.nx[0] || rindex[1] == grid.nx[1] || rindex[2] == grid.nx[2])
         sparse2(rindex,rindex,B,1.0,grid);
   }
   pb.epsilonp = 1.0;
   pb.epsilonm = 1.0;
   AMGsmall2(pb.psi,B,tmp.b,grid,2,0,grid.nx[0],grid.nx[0]*grid.nx[1]*grid.nx[2],
             1.0e-14,u,pb);
*/
   if (globlinsolve == 0)
   {
      cstart = clock ();
      BICGSTABsmall(pb.psi,tmp.A,tmp.b,grid,grid.nx[0]*grid.nx[1]*grid.nx[2],tollinsolve,
                    u,pb,tmp);
      cend = clock ();
      cout << "BICGSTAB clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
   }
   else if (globlinsolve == 1)
   {
      cstart = clock ();
      ILUsmall(tmp.M,tmp.A,grid,u,pb);
      cend = clock ();
      cout << "ILU clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
//      outputsparsesmall(beta2d,tmp.M,u,pb,grid);
//      preBICGSTABsmall(pb.psi,tmp.A,tmp.b,tmp.M,grid,grid.nx[0]*grid.nx[1]*grid.nx[2],
//                       tollinsolve,u,pb,tmp);
//      tmp.B = sparseelt2ptrmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
//      copyfromsmall(tmp.B,tmp.A,grid,u,pb);
      cstart = clock ();
      prerightBICGSTABsmall(pb.psi,tmp.A,tmp.b,tmp.M,grid,
                            grid.nx[0]*grid.nx[1]*grid.nx[2],tollinsolve,u,pb,tmp);
//      prerightBICGSTABsmall(pb.psi,tmp.B,tmp.b,tmp.M,grid,
//                            grid.nx[0]*grid.nx[1]*grid.nx[2],tollinsolve,u,pb,tmp);
      cend = clock ();
      cout << "preconditioned BICGSTAB clock time = " 
           << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
      gaussseidelsmall(pb.psi,tmp.A,tmp.b,grid,globGSsmooth,u,pb);
   }
   else if (globlinsolve == 2)
   {
    cstart = clock ();
     // AMGsmall2(pb.psi,tmp.A,tmp.b,grid,2,0,grid.nx[0],grid.nx[0]*grid.nx[1]*grid.nx[2],
     //           tollinsolve,u,pb);
      
      AMGsmall3(pb.psi,tmp.A,tmp.b,grid,2,0,grid.nx[0],grid.nx[0]*grid.nx[1]*grid.nx[2],
                tollinsolve,u,pb);
      cend = clock ();
      cout << "AMG clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
   }
   else if (globlinsolve == 3)
   {
      cstart = clock ();
      tmp.M = sparseelt2ptrmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
      ILUsmall(tmp.M,tmp.A,grid,u,pb);
      GMRESpreleftrestartsmall(pb.psi,tmp.A,tmp.b,tmp.M,grid,10,1.0e-12,u,pb,tmp);
//      GMRESrestartsmall(pb.psi,tmp.A,tmp.b,grid,20,1.0e-12,u,pb,tmp);
//      GMRESsmall(pb.psi,tmp.A,tmp.b,grid,10,tollinsolve,u,pb,tmp);
      cend = clock ();
      cout << "GMRES clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
   }else if(globlinsolve == 4 || globlinsolve == 5){
      cstart = clock ();
      amgsolve(pb.psi, tmp.A, tmp.b, grid,u, pb);
      cend = clock ();
      cout << "hypre clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
   }
   else
   {
      cout << "No linear solve option" << endl;
      exit(1);
   }
   if (globcheck)
   {
      // output(phi2d,pb.psi,grid.nx);
//      checkcim3Du(pb.psi,u,pb,grid);
      if(globcim !=7){
        checkanswer(pb.psi,u,pb,grid);
      }
      if (globcim == 5)
      {
         checkcim5Du(pb.psi,u,pb,grid);
         checkallDu(pb.psi,u,pb,grid);
      }
      else if (globcim == 345 || globcim == 0 ||globcim == 6)
      {
         if (globsmall)
            checkcim345Du(pb.psi,Dusmall,smallsize,u,pb,grid);
         else
            checkcim345Du(pb.psi,u,pb,grid);
         // checkcim345Dupm(pb.psi,u,pb,grid); //skip, not using Dusmall, slow
         // checkcim345jumpDu(pb.psi,u,pb,grid);
         // checkcim345DuAll(pb.psi,u,pb,grid);
      }
      // if (globcim == 5 || globcim == 345 || globcim == 0)
      // {
      //   // exit(1);
      // }
      // else{
      //   // getchar();
      // }
         
      else if (globcim == 2)
      {
        // checkDanswer2(pb.psi,u,pb,grid);
        checkcim12Du(pb.psi,u,pb,grid);
        if(globsmall){
          checkcim345Du(pb.psi,Dusmall,smallsize,u,pb,grid);
        }
      }
      else if (globcim == 4)
      {
        checkDanswer(pb.psi,u,pb,grid);
      }
      else if (globcim == 7)
      {
          double *** psi_true = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

          cout<<"Ghost state u error"<<endl;
          CheckErrGrid(pb.psi, sign_surf, pb, grid); // check with flipped surface          
          
          reconstruct(psi_true, u, pb.psi, sign_surf, grid);

          cout<<"Actual u error"<<endl;
          CheckErrGrid(psi_true, u, pb, grid); // check original surface 

          // Du err
          // cout<<"u error at interface"<<endl;
          // CheckIcimU(pb.psi, psi_true, sign_surf, u, pb, grid);

          cout<<"Du error at interface"<<endl;
          CheckIcimDu(pb.psi, psi_true, sign_surf, u, pb, grid);

          // CheckErrGeneral(psi_true, u, pb.psi, sign_surf, pb, grid);
          free_matrix(sign_surf, grid.nx[0],grid.nx[1],grid.nx[2]);
          free_matrix(psi_true, grid.nx[0],grid.nx[1],grid.nx[2] );
      }
      else
      {
        cerr<<"No recognized option in globcheck"<<endl;
        exit(1);
      }
      cout << "another ux err = " << globerrvec[0] << " " << globerrvec[1] << " " 
           << globerrvec[2] << endl;
      cout << "another uxx err = " << globerrvec2[0] << " " << globerrvec2[1] << " " 
           << globerrvec2[2] << endl;
      cout << "at " << globerrmx[0][0] << " " << globerrmx[0][1] << " " 
           << globerrmx[0][2] << endl;
      cout << "at " << globerrmx[1][0] << " " << globerrmx[1][1] << " " 
           << globerrmx[1][2] << endl;
      cout << "at " << globerrmx[2][0] << " " << globerrmx[2][1] << " " 
           << globerrmx[2][2] << endl;
      cout << "jump uxx err = " << globerrvec4[0] << " " << globerrvec4[1] << " " 
           << globerrvec4[2] << endl;
      cout << "exit globcheck in [getRHS]"<<endl;
      exit(1);
   }
//   exit(1);
// END CHANGE
   cstart = clock ();
   march.dorig = u;//original distance is u(surface)
   fmarch(march,-1.0,grid);
   cend = clock();
   cout << "fast march clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
//   output(phi2d,march.extend[0],grid.nx);
    if (globheatsmoothtime == 0.0){
      cout << "Do "<<globheatsmooth<<" steps of heatsmooth"<< endl;
      for (r = 1; r <= globheatsmooth; r++)
        advanceheat(march.extend[0],tmp,grid);
    }
   else
   {
      cout << "Do heatsmooth up to time "<<globheatsmoothtime<<endl;
      globheatsmoothtime = -grid.mindx*grid.mindx/(M_PI*M_PI)*log(grid.mindx);
      cout << "Doing " << (int) ceil(globheatsmoothtime/grid.mindx) << " heat steps" 
           << endl;
      advanceheat(march.extend[0],globheatsmoothtime,tmp,grid);
   }

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         for (s = 0, rindex[r] = max(tindex[r]-deg,0); s <= 2*deg; 
              s++, rindex[r] = min(max(tindex[r]-deg+s,0),grid.nx[r]))
            u1d[s] = evalarray(u,rindex);
         weno(dfp[r],dfn[r],u1d,deg,grid.dx[r],2*deg-1);
         rindex[r] = tindex[r];
      }
      if (fabs(evalarray(march.extend[0],tindex)) < grid.tol)
         setvalarray(rhs,tindex,0.0);
      else
      {
         if (evalarray(march.extend[0],tindex) >= grid.tol)
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
         setvalarray(rhs,tindex,-evalarray(march.extend[0],tindex)*grad);
      }
      if (maxvn < fabs(evalarray(march.extend[0],tindex))) //get maxvn
         maxvn = fabs(evalarray(march.extend[0],tindex));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

/*
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         Du[r] = 0.0;
         D2u[r][r] = -2.0*evalarray(u,tindex);
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = min(max(tindex[r]+s,0),grid.nx[r]);
            Du[r] += s*evalarray(u,rindex);
            D2u[r][r] += evalarray(u,rindex);
         }
         Du[r] /= 2.0*grid.dx[r];
         D2u[r][r] /= grid.dx[r]*grid.dx[r];
         rindex[r] = tindex[r];
      }

      for (r = 0; r < grid.dim; r++)
         for (s = r+1; s < grid.dim; s++)
         {
            D2u[r][s] = 0.0;
            for (j = -1; j <= 1; j += 2)
            {
               rindex[r] = min(max(tindex[r]+j,0),grid.nx[r]);
               for (k = -1; k <= 1; k += 2)
               {
                  rindex[s] = min(max(tindex[s]+k,0),grid.nx[s]);
                  D2u[r][s] += j*k*evalarray(u,rindex);
               }
            }
            D2u[r][s] /= 4.0*grid.dx[r]*grid.dx[s];
            rindex[r] = tindex[r];
            rindex[s] = tindex[s];
         }
      for (r = 0; r < grid.dim; r++)
         for (s = 0; s < r; s++)
            D2u[r][s] = D2u[s][r];

      value = 0.0;
      for (r = 0; r < grid.dim; r++)
         for (s = 0; s < grid.dim; s++)
            value -= Du[r]*D2u[r][s]*Du[s];
      grad2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         grad2 += Du[r]*Du[r];
      value /= max(grad2,tol);
      for (r = 0; r < grid.dim; r++)
         value += D2u[r][r];

      setvalarray(rhs,tindex,evalarray(rhs,tindex)+pb.gamma0*value);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
*/

// change time step here
//   grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim*pb.gamma0);
/*
   if (pb.gamma0 != 0.0 || maxvn >= 1.0)
      grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim*pb.gamma0+2.0*maxvn*grid.mindx);
   else
      grid.dt = grid.mindx/2.0;
*/

   if (globtestnum == 0)
      exactmaxvn = 3.0*sqrt(3.0)/4.0;
   else if (globtestnum == 1)
      exactmaxvn = fabs(2.0*(1.0-pb.epsilonp/pb.epsilonm)); //taken at r = 1, 429.48
   else if (globtestnum == 2)
      exactmaxvn = 2.0*exp(1.0)*fabs(1.0-pb.epsilonp/pb.epsilonm); //taken at r = 1, 158
   else if (globtestnum == 3 || globtestnum == 4) 
//      exactmaxvn = fabs(1.0-pb.epsilonp/pb.epsilonm);
      exactmaxvn = fabs(1.0-pb.epsilonp/pb.epsilonm)/
                   sqrt(1.0+fabs(1.0-pb.epsilonp/pb.epsilonm));
   cout << "calculated maxvn = " << maxvn << " while exact is " << exactmaxvn << endl;

   // if (maxvn > exactmaxvn)
   // {
   //    cout << "Exact maxvn is not max" << endl;
   //    exit(1);
   // }
   // else
   //    maxvn = exactmaxvn;

   if (!globdtdx2)
      grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim*pb.gamma0+2.0*maxvn*grid.mindx);
   else if (globtestnum == 0);
   else if (globtestnum == 1)
      grid.dt = grid.mindx*grid.mindx/22.0;
   else if (globtestnum == 2)
      grid.dt = grid.mindx*grid.mindx/58.0;
   else if (globtestnum == 3 || globtestnum == 4)
      grid.dt = grid.mindx*grid.mindx/1.2; //extra 1/2 for investigation

   if (grid.t+grid.dt > grid.tfinal)
      grid.dt = grid.tfinal-grid.t;
   cout << "   dt = " << grid.dt << " because max vn = " << maxvn << endl;
  checkwithexactvn(march.extend[0], u, pb, grid); // check vn with exact
  #if 1//no Dusmall in linearsystem4, pointer being freed was not allocated 
   SparseElt *current, *current2;
   for (i = 0; i < smallsize; i++)
   {
      for (current = Dusmall[i].head, current2 = current; current != NULL; 
           current = (*current).next, delete current2, current2 = current);
      delete [] Dusmall[i].info;
   }
   delete [] Dusmall;
   #endif
   
}
// with a
void getRHS(double ***rhs, double ***u, double ***a, PBData &pb, MarchStruct &march, 
            TempStruct &tmp, GridData &grid)
{
   double dfp[grid.dim], dfn[grid.dim], df[grid.dim];
   double grad, maxvn = 0.0, exactmaxvn;
   int deg = 3;
   double u1d[2*deg+1];
   int i, r, s, j, k;
   int tindex[grid.dim], rindex[grid.dim];
   double Du[grid.dim], D2u[grid.dim][grid.dim];
   double value, grad2, tol = grid.tol;
   // ofstream phi2d("phi2d.dat", ios::out);
   // ofstream beta2d("beta2d.dat", ios::out);
   clock_t cstart, cend;

   for (r = 0; r < grid.dim; r++)
   {
      globerrvec[r] = 0.0;
      globerrvec2[r] = 0.0;
      globerrvec3[r] = 0.0;
      globerrvec4[r] = 0.0;
      for (s = 0; s < grid.dim; s++)
         globerrmx[r][s] = 0.0;
   }

   cout << "in getRHS" << endl;
// CHANGE HERE
// solves the Poisson-Boltzman equation and puts the answer in the 3D array pb.psi
   cstart = clock ();
   StorageStruct *Dusmall;
   int smallsize = 0;
   linearsystem6(tmp.A,tmp.b,Dusmall,smallsize,a,u,pb,grid);
   march.Dusmall = Dusmall;
   march.smallsize = smallsize;
   cend = clock ();
   cout << "linear system clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC 
        << endl;
   cout << "equation error = " << globerr << " at " << globerrvec3[0] << " " 
        << globerrvec3[1] << " " << globerrvec3[2] << endl;
   if (globlinsolve == 0)
   {
      cstart = clock ();
      BICGSTABsmall(pb.psi,tmp.A,tmp.b,grid,grid.nx[0]*grid.nx[1]*grid.nx[2],tollinsolve,
                    a,u,pb,tmp);
      cend = clock ();
      cout << "BICGSTAB clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
   }
   else if (globlinsolve == 1)
   {
      cstart = clock ();
      ILUsmall(tmp.M,tmp.A,grid,a,u,pb);
      cend = clock ();
      cout << "ILU clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
      cstart = clock ();
      prerightBICGSTABsmall(pb.psi,tmp.A,tmp.b,tmp.M,grid,
                            grid.nx[0]*grid.nx[1]*grid.nx[2],tollinsolve,a,u,pb,tmp);
      cend = clock ();
      cout << "preconditioned BICGSTAB clock time = " 
           << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
      gaussseidelsmall(pb.psi,tmp.A,tmp.b,grid,globGSsmooth,a,u,pb);
   }
   else if (globlinsolve == 2)
   {
      cstart = clock ();
      AMGsmall3(pb.psi,tmp.A,tmp.b,grid,2,0,grid.nx[0],grid.nx[0]*grid.nx[1]*grid.nx[2],
                tollinsolve,a,u,pb);
      cend = clock ();
      cout << "AMG clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
   }
   else
   {
      cout << "No linear solve option" << endl;
      exit(1);
   }
   if (globcheck)
   {
      // output(phi2d,pb.psi,grid.nx);
      checkanswer(pb.psi,u,pb,grid);
      if (globcim == 5)
      {
         checkcim5Du(pb.psi,u,pb,grid);
         checkallDu(pb.psi,u,pb,grid);
      }
      else if (globcim == 345 || globcim == 0 || globcim == 6)
      {
         if (globsmall)
            checkcim345Du(pb.psi,Dusmall,smallsize,u,pb,grid);
         else
            checkcim345Du(pb.psi,u,pb,grid);
      }
      // if (globcim == 5 || globcim == 345 || globcim == 0)
      //    {exit(1);}
      // else{
      //   // getchar();
      // }
         
      if (globcim == 2)
         checkDanswer2(pb.psi,u,pb,grid);
      else if (globcim == 4)
         checkDanswer(pb.psi,u,pb,grid);
      cout << "another ux err = " << globerrvec[0] << " " << globerrvec[1] << " " 
           << globerrvec[2] << endl;
      cout << "another uxx err = " << globerrvec2[0] << " " << globerrvec2[1] << " " 
           << globerrvec2[2] << endl;
      cout << "at " << globerrmx[0][0] << " " << globerrmx[0][1] << " " 
           << globerrmx[0][2] << endl;
      cout << "at " << globerrmx[1][0] << " " << globerrmx[1][1] << " " 
           << globerrmx[1][2] << endl;
      cout << "at " << globerrmx[2][0] << " " << globerrmx[2][1] << " " 
           << globerrmx[2][2] << endl;
      cout << "jump uxx err = " << globerrvec4[0] << " " << globerrvec4[1] << " " 
           << globerrvec4[2] << endl;
      cout << "exit globcheck in [getRHS]"<<endl;
      exit(1);
   }
//   exit(1);
// END CHANGE
   cstart = clock ();
   march.dorig = u;
   fmarch(march,-1.0,grid);
   cend = clock();
   cout << "fast march clock time = " << (double) (cend-cstart)/CLOCKS_PER_SEC << endl;
//   output(phi2d,march.extend[0],grid.nx);
   if (globheatsmoothtime == 0.0)
      for (r = 1; r <= globheatsmooth; r++)
         advanceheat(march.extend[0],tmp,grid);
   else
   {
      globheatsmoothtime = -grid.mindx*grid.mindx/(M_PI*M_PI)*log(grid.mindx);
      cout << "Doing " << (int) ceil(globheatsmoothtime/grid.mindx) << " heat steps" 
           << endl;
      advanceheat(march.extend[0],globheatsmoothtime,tmp,grid);
   }

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         for (s = 0, rindex[r] = max(tindex[r]-deg,0); s <= 2*deg; 
              s++, rindex[r] = min(max(tindex[r]-deg+s,0),grid.nx[r]))
            u1d[s] = evalarray(u,rindex);
         weno(dfp[r],dfn[r],u1d,deg,grid.dx[r],2*deg-1);
         rindex[r] = tindex[r];
      }
      if (fabs(evalarray(march.extend[0],tindex)) < grid.tol)
         setvalarray(rhs,tindex,0.0);
      else
      {
         if (evalarray(march.extend[0],tindex) >= grid.tol)
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
         setvalarray(rhs,tindex,-evalarray(march.extend[0],tindex)*grad);
      }
      if (maxvn < fabs(evalarray(march.extend[0],tindex)))
         maxvn = fabs(evalarray(march.extend[0],tindex));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

/*
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         Du[r] = 0.0;
         D2u[r][r] = -2.0*evalarray(u,tindex);
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = min(max(tindex[r]+s,0),grid.nx[r]);
            Du[r] += s*evalarray(u,rindex);
            D2u[r][r] += evalarray(u,rindex);
         }
         Du[r] /= 2.0*grid.dx[r];
         D2u[r][r] /= grid.dx[r]*grid.dx[r];
         rindex[r] = tindex[r];
      }

      for (r = 0; r < grid.dim; r++)
         for (s = r+1; s < grid.dim; s++)
         {
            D2u[r][s] = 0.0;
            for (j = -1; j <= 1; j += 2)
            {
               rindex[r] = min(max(tindex[r]+j,0),grid.nx[r]);
               for (k = -1; k <= 1; k += 2)
               {
                  rindex[s] = min(max(tindex[s]+k,0),grid.nx[s]);
                  D2u[r][s] += j*k*evalarray(u,rindex);
               }
            }
            D2u[r][s] /= 4.0*grid.dx[r]*grid.dx[s];
            rindex[r] = tindex[r];
            rindex[s] = tindex[s];
         }
      for (r = 0; r < grid.dim; r++)
         for (s = 0; s < r; s++)
            D2u[r][s] = D2u[s][r];

      value = 0.0;
      for (r = 0; r < grid.dim; r++)
         for (s = 0; s < grid.dim; s++)
            value -= Du[r]*D2u[r][s]*Du[s];
      grad2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         grad2 += Du[r]*Du[r];
      value /= max(grad2,tol);
      for (r = 0; r < grid.dim; r++)
         value += D2u[r][r];

      setvalarray(rhs,tindex,evalarray(rhs,tindex)+pb.gamma0*value);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
*/

// change time step here
//   grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim*pb.gamma0);
/*
   if (pb.gamma0 != 0.0 || maxvn >= 1.0)
      grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim*pb.gamma0+2.0*maxvn*grid.mindx);
   else
      grid.dt = grid.mindx/2.0;
*/

   if (globtestnum == 0)
      exactmaxvn = 3.0*sqrt(3.0)/4.0;
   else if (globtestnum == 1)
      exactmaxvn = fabs(2.0*(1.0-pb.epsilonp/pb.epsilonm));
   else if (globtestnum == 2)
      exactmaxvn = 2.0*exp(1.0)*fabs(1.0-pb.epsilonp/pb.epsilonm);
   else if (globtestnum == 3 || globtestnum == 4)
//      exactmaxvn = fabs(1.0-pb.epsilonp/pb.epsilonm);
      exactmaxvn = fabs(1.0-pb.epsilonp/pb.epsilonm)/
                   sqrt(1.0+fabs(1.0-pb.epsilonp/pb.epsilonm));
   cout << "calculated maxvn = " << maxvn << " while exact is " << exactmaxvn << endl;

   // if (maxvn > exactmaxvn)// do not stop
   // {
   //    cout << "Exact maxvn is not max" << endl;
   //    exit(1);
   // }
   // else
   //    maxvn = exactmaxvn;

   if (!globdtdx2)
      grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim*pb.gamma0+2.0*maxvn*grid.mindx);
   else if (globtestnum == 0);
   else if (globtestnum == 1)
      grid.dt = grid.mindx*grid.mindx/22.0;
   else if (globtestnum == 2)
      grid.dt = grid.mindx*grid.mindx/58.0;
   else if (globtestnum == 3 || globtestnum == 4)
      grid.dt = grid.mindx*grid.mindx/1.2;

   if (grid.t+grid.dt > grid.tfinal)
      grid.dt = grid.tfinal-grid.t;
   cout << "   dt = " << grid.dt << " because max vn = " << maxvn << endl;
   checkwithexactvn(march.extend[0], u, pb, grid); // check vn with exact
   #if 1
   SparseElt *current, *current2;
   for (i = 0; i < smallsize; i++)
   {
      for (current = Dusmall[i].head, current2 = current; current != NULL; 
           current = (*current).next, delete current2, current2 = current);
      delete [] Dusmall[i].info;
   }
   delete [] Dusmall;
   #endif
}
// globexactmotion, thevn is calculated from location of interface, not extended by fast marching
void getRHShyp(double ***rhs, double ***u, PBData &pb, MarchStruct &march, 
               TempStruct &tmp, GridData &grid)
{
   double dfp[grid.dim], dfn[grid.dim], df[grid.dim];
   double grad, maxvn = 0.0, thevn, theradius, x[grid.dim];
   int deg = 3;
   double u1d[2*deg+1];
   int i, r, s, j, k;
   int tindex[grid.dim], rindex[grid.dim];
   double Du[grid.dim], D2u[grid.dim][grid.dim];
   double value, grad2, tol = grid.tol;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      sub2coord(x,tindex,grid);
      theradius = 0;
      for (r = 0; r < grid.dim; r++)
         theradius += x[r]*x[r];
      theradius = sqrt(theradius);
      thevn = 4.0*theradius/((1.0+theradius*theradius)*(1.0+theradius*theradius));

      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         for (s = 0, rindex[r] = max(tindex[r]-deg,0); s <= 2*deg; 
              s++, rindex[r] = min(max(tindex[r]-deg+s,0),grid.nx[r]))
            u1d[s] = evalarray(u,rindex);
         weno(dfp[r],dfn[r],u1d,deg,grid.dx[r],2*deg-1);
         rindex[r] = tindex[r];
      }
      if (fabs(thevn) < grid.tol)
         setvalarray(rhs,tindex,0.0);
      else
      {
         if (thevn > 0.0)
            for (r = 0; r < grid.dim; r++)
            {
               dfp[r] = (dfp[r] < 0.0) ? -dfp[r] : 0.0; // max(-dfp,0)
               dfn[r] = (dfn[r] > 0.0 ) ? dfn[r] : 0.0; // max( dfn,0)
               df[r]  = (dfp[r] > dfn[r]) ? dfp[r] : dfn[r]; // max( max(-dfp,0), max( dfn,0) ) = max(dfn, -dfp, 0)
            }
         else
            for (r = 0; r < grid.dim; r++)
            {
               dfp[r] = (dfp[r] > 0.0) ? dfp[r] : 0.0;
               dfn[r] = (dfn[r] < 0.0) ? -dfn[r] : 0.0;
               df[r]  = (dfp[r] > dfn[r]) ? dfp[r] : dfn[r]; // max(dfp, -dfn, 0)
            }
         grad = 0.0;
         for (r = 0; r < grid.dim; r++)
            grad += df[r]*df[r];
         grad = sqrt(grad);
         setvalarray(rhs,tindex,-thevn*grad);
      }
      if (maxvn < fabs(thevn))
         maxvn = fabs(thevn);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   maxvn = 3.0*sqrt(3.0)/4.0;
   grid.dt = grid.mindx/(2.0*maxvn);
   if (grid.t+grid.dt > grid.tfinal)
      grid.dt = grid.tfinal-grid.t;
   cout << "   dt = " << grid.dt << " because max vn = " << maxvn << endl;
}
//globexactvn, u is surface, 
void getRHSexactvn(double ***rhs, double ***u, PBData &pb, MarchStruct &march, 
                   TempStruct &tmp, GridData &grid)
{
   double dfp[grid.dim], dfn[grid.dim], df[grid.dim];
   double grad, maxvn = 0.0, thevn, theradius, x[grid.dim];
   int deg = 3;
   double u1d[2*deg+1];
   int i, r, s, j, k;
   int tindex[grid.dim], rindex[grid.dim];
   double Du[grid.dim], D2u[grid.dim][grid.dim];
   double value, grad2, tol = grid.tol;

   march.dorig = u;
   fmarch(march,-1.0,grid);
   if (globtestnum == 0)
      maxvn = 3.0*sqrt(3.0)/4.0;
   // maxvn = 2.0*fabs(2.0*(1.0-pb.epsilonp/pb.epsilonm));
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
/*
      sub2coord(x,tindex,grid);
      theradius = 0.0;
      for (r = 0; r < grid.dim; r++)
         theradius += x[r]*x[r];
      theradius = sqrt(theradius);
      thevn = 2.0*(1.0-pb.epsilonp/pb.epsilonm)*theradius;
*/
      thevn = evalarray(march.extend[0],tindex);
      if (fabs(thevn) > maxvn)
      {
         if (thevn > 0.0)
            thevn = maxvn;
         else
            thevn = -maxvn;
      }

      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         for (s = 0, rindex[r] = max(tindex[r]-deg,0); s <= 2*deg; 
              s++, rindex[r] = min(max(tindex[r]-deg+s,0),grid.nx[r]))
            u1d[s] = evalarray(u,rindex);
         weno(dfp[r],dfn[r],u1d,deg,grid.dx[r],2*deg-1);
         rindex[r] = tindex[r];
      }
      if (fabs(thevn) < grid.tol)
         setvalarray(rhs,tindex,0.0);
      else
      {
         if (thevn > 0.0)
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
         setvalarray(rhs,tindex,-thevn*grad);
      }
      if (maxvn < fabs(thevn))
         maxvn = fabs(thevn);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   grid.dt = grid.mindx/(2.0*maxvn);
   if (grid.t+grid.dt > grid.tfinal)
      grid.dt = grid.tfinal-grid.t;
   cout << "   dt = " << grid.dt << " because max vn = " << maxvn << endl;
}

// motion by mean cuvature
void getRHScurv(double ***rhs, double ***u, TempStruct &tmp, GridData &grid)
{
   int i, r, s, j, k;
   int tindex[grid.dim], rindex[grid.dim];
   double Du[grid.dim], D2u[grid.dim][grid.dim];
   double value, grad2, tol = grid.tol;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         Du[r] = 0.0;
         D2u[r][r] = -2.0*evalarray(u,tindex);
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = min(max(tindex[r]+s,0),grid.nx[r]);
            Du[r] += s*evalarray(u,rindex);
            D2u[r][r] += evalarray(u,rindex);
         }
         Du[r] /= 2.0*grid.dx[r];
         D2u[r][r] /= grid.dx[r]*grid.dx[r];
         rindex[r] = tindex[r];
      }

      for (r = 0; r < grid.dim; r++)
         for (s = r+1; s < grid.dim; s++)
         {
            D2u[r][s] = 0.0;
            for (j = -1; j <= 1; j += 2)
            {
               rindex[r] = min(max(tindex[r]+j,0),grid.nx[r]);
               for (k = -1; k <= 1; k += 2)
               {
                  rindex[s] = min(max(tindex[s]+k,0),grid.nx[s]);
                  D2u[r][s] += j*k*evalarray(u,rindex);
               }
            }
            D2u[r][s] /= 4.0*grid.dx[r]*grid.dx[s];
            rindex[r] = tindex[r];
            rindex[s] = tindex[s];
         }
      for (r = 0; r < grid.dim; r++)
         for (s = 0; s < r; s++)
            D2u[r][s] = D2u[s][r];

      value = 0.0;
      for (r = 0; r < grid.dim; r++)
         for (s = 0; s < grid.dim; s++)
            value -= Du[r]*D2u[r][s]*Du[s];
      grad2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         grad2 += Du[r]*Du[r];
      value /= max(grad2,tol);
      for (r = 0; r < grid.dim; r++)
         value += D2u[r][r];

      setvalarray(rhs,tindex,value);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim);
}




double getpsivac(double *x, PBData &pb)
{
   double radius, value = 0.0;
   int i, j;

   for (i = 0; i < pb.N; i++)
   {
      radius = 0.0;
      for (j = 0; j < pb.dim; j++)
         radius += (x[j]-pb.x[i][j])*(x[j]-pb.x[i][j]);
      radius = sqrt(radius);
      value += pb.Q[i]/(4.0*M_PI*pb.epsilonm*radius);
   }

   return value;
}

void getDpsivac(double *Dpsi, double *x, PBData &pb)
{
   double radius;
   int i, j;

   for (j = 0; j < pb.dim; j++)
      Dpsi[j] = 0.0;

   for (i = 0; i < pb.N; i++)
   {
      radius = 0.0;
      for (j = 0; j < pb.dim; j++)
         radius += (x[j]-pb.x[i][j])*(x[j]-pb.x[i][j]);
      for (j = 0; j < pb.dim; j++)
         Dpsi[j] += -pb.Q[i]/(4.0*M_PI*pb.epsilonm*radius*sqrt(radius))*
                     (x[j]-pb.x[i][j]);
   }
}

double getDpsivacn(int *index, double ***S, PBData &pb, GridData &grid)
{
   double value, x[grid.dim], Dpsi[grid.dim], normal[grid.dim];
   int r, s, tindex[grid.dim], rindex[grid.dim];

   sub2coord(x,index,grid);
   getDpsivac(Dpsi,x,pb);
   
   for (r = 0; r < grid.dim; r++)
      tindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
   {
      normal[r] = 0.0;
      for (s = -1; s <= 1; s += 2)
      {
         tindex[r] = index[r]+s;
         normal[r] += s*evalarray(S,tindex);
      }
      normal[r] /= 2.0*grid.dx[r];
      tindex[r] = index[r];
   }

   value = 0.0;
   for (r = 0; r < grid.dim; r++)
      value += normal[r]*normal[r];
   value = sqrt(value);
   for (r = 0; r < grid.dim; r++)
      normal[r] /= value;

   value = 0.0;
   for (r = 0; r < grid.dim; r++)
      value += Dpsi[r]*normal[r];

   return value;
}

void getD2psivac(double **D2psi, double *x, PBData &pb)
{
   double radius;
   int i, j, m;

   for (j = 0; j < pb.dim; j++)
      for (m = 0; m < pb.dim; m++)
         D2psi[j][m] = 0.0;

   for (i = 0; i < pb.N; i++)
   {
      radius = 0.0;
      for (j = 0; j < pb.dim; j++)
         radius += (x[j]-pb.x[i][j])*(x[j]-pb.x[i][j]);
      for (j = 0; j < pb.dim; j++)
      {
         for (m = 0; m < pb.dim; m++)
            D2psi[j][m] += 3.0*pb.Q[i]/(4.0*M_PI*pb.epsilonm*radius*radius*
                                        sqrt(radius))*
                           (x[j]-pb.x[i][j])*(x[m]-pb.x[i][m]);
         D2psi[j][j] += -pb.Q[i]/(4.0*M_PI*pb.epsilonm*radius*sqrt(radius));
      }
   }
}

double Bval(double s, PBData &pb)
{
   double value = 0.0;
   int i;

   for (i = 0; i < pb.N; i++)
      value += pb.c[i]*(exp(-pb.beta*pb.epsilone[i]*s)-1);
   value /= pb.beta;

   return value;
}

double Bprime(double s, PBData &pb)
{
   double value = 0.0;
   int i;

   for (i = 0; i < pb.N; i++)
      value -= pb.epsilone[i]*pb.c[i]*exp(-pb.beta*pb.epsilone[i]*s);

   return value;
}

double B2prime(double s, PBData &pb)
{
   double value = 0.0;
   int i;

   for (i = 0; i < pb.N; i++)
      value += pb.c[i]*pb.beta*pb.epsilone[i]*pb.epsilone[i]*
               exp(-pb.beta*pb.epsilone[i]*s);

   return value;
}

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
// get maximum error of u at grid points
void checkanswer(double ***u, double ***S, PBData &pb, GridData &grid)
{
   int i, s, tindex[grid.dim], rindex[grid.dim];
   double theerr = 0.0, tmperr;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (evalarray(S,tindex) < 0.0)
         tmperr = evalarray(u,tindex)-getu(tindex,0,0,0.0,-1,grid);
      else
         tmperr = evalarray(u,tindex)-getu(tindex,0,0,0.0,1,grid);
      if (fabs(tmperr) > fabs(theerr))
      {
         theerr = tmperr;
         for (s = 0; s < grid.dim; s++)
            rindex[s] = tindex[s];
      }

      if (globwriteerr){
        outfile_uerr<<tindex[0]<<","<<tindex[1]<<","<<tindex[2]<<",";
        outfile_uerr <<setprecision(12)<<tmperr<< endl;;
      }


      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
    
  cout<<"[checkanswer]"<<endl;
   cout << "Error is " << theerr << " at " << rindex[0] << " " << rindex[1] << " "  
        << rindex[2] << endl;
}

void checkDanswer(double ***u, double ***S, PBData &pb, GridData &grid)
{
   int gamma[grid.dim][2];
   int i, r, s, j, k, thestatus;
   int tindex[grid.dim], rindex[grid.dim];
   double uxerr[grid.dim];
   char except;
   int tempsk2[4];

   for (r = 0; r < grid.dim; r++)
      uxerr[r] = 0.0;

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
      thestatus = getstatus3(S,tindex,grid);

      if (thestatus == 1);
      else if (thestatus == 2 || thestatus == 3)
      {
         except = 0;
         for (r = 0; r < grid.dim && !except; r++)
            for (s = r+1; s < grid.dim && !except; s++)
               if (!yessk2(tempsk2,r,s,tindex,S,grid))
                  except = 1;
         if (!except)
            checkDanswer(uxerr,tindex,gamma,u,S,pb,grid);
         else;
/*
         if (thestatus == 2)
            checkDanswer(uxerr,tindex,gamma,u,S,pb,grid);
*/
      }
      else;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "[checkDanswer]DError is " << uxerr[0] << " " << uxerr[1] << " " << uxerr[2] << endl;
}

void checkDanswer(double *uxerr, int *index, int gamma[][2], double ***u, double ***S, 
                  PBData &pb, GridData &grid)
{
   int r, s, t, i, m, n;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];
   double d0[grid.dim];
   double ***dcoef[grid.dim];
//   double ***D1[grid.dim], ***D2[grid.dim][grid.dim];
   double ***D1[grid.dim], ***D2[grid.dim][3];
   double **LU, **G, value;
//   int sk2[grid.dim][grid.dim][4];
   int sk2[grid.dim][3][4];
   int PLR[grid.dim], PLC[grid.dim];
   double alpha[grid.dim], beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[grid.dim];
   int sk[grid.dim];
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;
   double jumpD1u, ***jumpD1ucoef, jumpD1uxxcoef[grid.dim];
   double jumpD2u, ***jumpD2ucoef, jumpD2uxxcoef[grid.dim];
   double ****D1ucoef[grid.dim], ***D1uxxcoef;
   double uxx[grid.dim], tmperr[grid.dim], realuxx[grid.dim];

   double zerouxxcoef[grid.dim];
   double ***zeroucoef;
   zeroucoef = matrix(2,2,2);
   for (s = 0; s < grid.dim; s++) 
      zerouxxcoef[s] = 0.0;
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      setvalarray(zeroucoef,sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;
   value = 0.0;

   LU = matrix(grid.dim-1,grid.dim-1);
   G = matrix(grid.dim-1,grid.dim-1);
   jumpD1ucoef = matrix(2,2,2);
   jumpD2ucoef = matrix(2,2,2);
   for (r = 0; r < grid.dim; r++)
   {
      D1[r] = matrix(2,2,2);
      for (s = 0; s < grid.dim; s++)
         D2[r][s] = matrix(2,2,2);
      dcoef[r] = matrix(2,2,2);
   }
   for (r = 0; r < grid.dim; r++)
      D1ucoef[r] = matrix(grid.dim-1,2,2,2);
   D1uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);

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
      sk[r] = gamma[r][1]-gamma[r][0];
// get sk2
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
         getsk2(sk2[m][n],m,n,index,S,grid);

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;
   for (r = 0; r < grid.dim; r++)
      if (abs(sk[r]) == 1)
      {
         rindex[r] = index[r]+sk[r];
// getting first and second derivatives
         getD1(D1,sk,grid);
         getD2(D2,sk2,grid);
// getting interface info
         getinterfaceinfo(alpha[r],tangent,normal,S,index,rindex,grid);
         beta = 1.0-alpha[r];
// get Du 2nd order
         getDu(D1ucoef[r],D1uxxcoef[r],index,r,sk[r],alpha[r],thesign,sk,D1,D2,grid);
// getting jump of uxx in r
         gettau(tau,index,r,sk[r],alpha[r],grid);
         getjumpux(jumpD1u,jumpD1ucoef,jumpD1uxxcoef,index,r,sk,alpha[r],thesign,
                   normal,tangent,D1ucoef[r],D1uxxcoef[r],S,pb,grid);
         if (globorder == 1)
            getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,index,r,sk,alpha[r],thesign,
                       normal,D1,D2,S,pb,grid);
         else
            getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,index,r,sk,alpha[r],thesign,
                       normal,D1ucoef[r],D1uxxcoef[r],D2,S,pb,grid);
         if (fabs(evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,index,S,grid)-
                  (getD2u(index,r,r,r,sk[r],alpha[r],1,grid)-
                   getD2u(index,r,r,r,sk[r],alpha[r],-1,grid))) > globerrvec4[r])
            globerrvec4[r] = fabs(evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,
                                           index,S,grid)-
                                  (getD2u(index,r,r,r,sk[r],alpha[r],1,grid)-
                                   getD2u(index,r,r,r,sk[r],alpha[r],-1,grid)));
// form d0 and dcoef's rth entry 
         d0[r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+sk[r]*beta*jumpD1u/grid.dx[r]+
                          0.5*beta*beta*jumpD2u);

         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            setvalarray(dcoef[r],sindex,
                        thesign*(sk[r]*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                 0.5*beta*beta*evalarray(jumpD2ucoef,sindex))-
                        sk[r]*(alpha[r]+beta)*evalarray(D1ucoef[r][r],sindex)/
                        grid.dx[r]);

            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;
         setvalarray(dcoef[r],sindex,
                     evalarray(dcoef[r],sindex)-1.0/(grid.dx[r]*grid.dx[r]));
         sindex[r] = 1+sk[r];
         setvalarray(dcoef[r],sindex,
                     evalarray(dcoef[r],sindex)+1.0/(grid.dx[r]*grid.dx[r]));
         sindex[r] = 1;
// form matrix's rth row
         for (m = 0; m < grid.dim; m++)
            G[r][m] = -thesign*(sk[r]*beta*jumpD1uxxcoef[m]/grid.dx[r]+
                                0.5*beta*beta*jumpD2uxxcoef[m])+
                      sk[r]*(alpha[r]+beta)*D1uxxcoef[r][r][m]/grid.dx[r];
         G[r][r] += 0.5*(beta*beta-alpha[r]*alpha[r]);

         rindex[r] = index[r];
      }
      else
      {
// if interior point
         d0[r] = 0.0;
         for (s = 0; s < grid.dim; s++)
            G[r][s] = 0.0;
         G[r][r] = 1.0;
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            setvalarray(dcoef[r],sindex,0.0);

            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;
         setvalarray(dcoef[r],sindex,-2.0/(grid.dx[r]*grid.dx[r]));
         for (s = -1; s <= 1; s += 2)
         {
            sindex[r] = 1+s;
            setvalarray(dcoef[r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
         }
         sindex[r] = 1;
      }

   gecp0(LU,PLR,PLC,G,grid.dim-1,grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,grid.dim-1);
   for (s = 0; s < grid.dim; s++)
   {
      uxx[s] = temp[s];
      realuxx[s] = temp[s];
   }

   for (s = 0; s < grid.dim; s++)
      tindex[s] = index[s];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] < 3)
   {
      for (n = 0; n < grid.dim; n++)
         temp[n] = evalarray(dcoef[n],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-1+sindex[s];
      for (s = 0; s < grid.dim; s++)
      {
         uxx[s] += temp[s]*evalarray(u,tindex);
         if (evalarray(S,tindex) < 0.0)
            realuxx[s] += temp[s]*getu(tindex,0,0,0.0,-1,grid);
         else
            realuxx[s] += temp[s]*getu(tindex,0,0,0.0,1,grid);
      }

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 1;

   double anoDu[grid.dim];
   for (r = 0; r < grid.dim; r++)
      if (abs(sk[r]) == 1)
      {
         for (s = 0; s < grid.dim; s++)
         {
            tmperr[s] = 0.0;
            for (m = 0; m < grid.dim; m++)
               tmperr[s] += D1uxxcoef[r][s][m]*uxx[m];
         }

         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            for (s = 0; s < grid.dim; s++)
               tindex[s] = index[s]-1+sindex[s];
            for (s = 0; s < grid.dim; s++)
               tmperr[s] += evalarray(D1ucoef[r][s],sindex)*evalarray(u,tindex);

            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;

         for (s = 0; s < grid.dim; s++)
         {
            tmperr[s] -= getDu(index,s,r,sk[r],alpha[r],thesign,grid);
            if (fabs(tmperr[s]) > fabs(uxerr[s]))
               uxerr[s] = tmperr[s];
         }

         for (s = 0; s < grid.dim; s++)
         {
            tmperr[s] = evalcoef(0.0,D1ucoef[r][s],D1uxxcoef[r][s],index,S,grid)-
                        getDu(index,s,r,sk[r],alpha[r],thesign,grid);
            if (fabs(tmperr[s]) > fabs(globerrvec[s]))
               globerrvec[s] = tmperr[s];
            tmperr[s] = uxx[s]-getD2u(index,s,s,r,sk[r],alpha[r],thesign,grid);
            if (fabs(tmperr[s]) > fabs(globerrvec2[s]))
            {
               globerrvec2[s] = tmperr[s];
               for (m = 0; m < grid.dim; m++)
                  globerrmx[s][m] = index[m];
            }
         }
      }

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]) 
   {
      cout << "IN CHECK ERROR" << endl;
      cout << uxx[0] << " " << uxx[1] << " " << uxx[2] << endl;
      cout << realuxx[0] << " " << realuxx[1] << " " << realuxx[2] << endl;
      cout << getD2u(index,0,0,0,0,0.0,thesign,grid) << " "
           << getD2u(index,1,1,0,0,0.0,thesign,grid) << " "
           << getD2u(index,2,2,0,0,0.0,thesign,grid) << endl;
   }

   free_matrix(zeroucoef,2,2,2);
   free_matrix(LU,grid.dim-1,grid.dim-1);
   free_matrix(G,grid.dim-1,grid.dim-1);
   free_matrix(jumpD1ucoef,2,2,2);
   free_matrix(jumpD2ucoef,2,2,2);
   for (r = 0; r < grid.dim; r++)
   {
      free_matrix(D1[r],2,2,2);
      for (s = 0; s < grid.dim; s++)
         free_matrix(D2[r][s],2,2,2);
      free_matrix(dcoef[r],2,2,2);
   }
   for (r = 0; r < grid.dim; r++)
      free_matrix(D1ucoef[r],grid.dim-1,2,2,2);
   free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
}
//check Du from cim2 solution
void checkDanswer2(double ***u, double ***S, PBData &pb, GridData &grid)
{
   int gamma[grid.dim][2];
   int i, r, s, j, k, thestatus;
   int tindex[grid.dim], rindex[grid.dim];
   double uxerr[grid.dim];

   for (r = 0; r < grid.dim; r++)
      uxerr[r] = 0.0;

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
      thestatus = getstatus3(S,tindex,grid);

      if (thestatus == 1);//interior
      else if (thestatus == 2)//can use CIM2
         checkDanswer2(uxerr,tindex,gamma,u,S,pb,grid);
      else if (thestatus == 3)
        cout<<"cim1 check not implemented"<<endl;//exceptional
      else
        cout<<"Undefined status in checkDanswer2"<<endl;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "[checkDanswer2]DError is " << uxerr[0] << " " << uxerr[1] << " " << uxerr[2] << endl;
}
// check cim2 error, store in glob err vec
void checkDanswer2(double *uxerr, int *index, int gamma[][2], double ***u, double ***S, 
                   PBData &pb, GridData &grid)
{
   int r, s, i, j, k, tmps;
   int rindex[grid.dim], tindex[grid.dim];
   double ****f; 
   double **LU, **M, value, Sval; 
   int PLR[grid.dim], PLC[grid.dim];
   double alpha[grid.dim], beta, tangent[grid.dim], normal[grid.dim];
   double bk, rthere, rhere, ethere, ehere, ehat, C;
   double temp[grid.dim], a[4];
   int sk[grid.dim];
// ADDING
   double sigma, d[grid.dim], thesign, tau, Dtau[grid.dim];

   LU = matrix(grid.dim-1,grid.dim-1); //3 by 3 matrix
   M = matrix(grid.dim-1,grid.dim-1); //3 by 3 matrix
   f = matrix(grid.dim-1,4,4,4); //3x5x5x5 matrix : coefficient of u for DDu in dim r. current point is [2,2,2].

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
      sk[r] = gamma[r][1]-gamma[r][0]; //sk[r] = s = +-1, interface in [tindex[r],tindex[r+s]]

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

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      tindex[r] = 2;
   for (r = 0; r < grid.dim; r++)
   {
      rindex[r] = index[r]+sk[r];
      getinterfaceinfo(alpha[r],tangent,normal,S,index,rindex,grid);
      rindex[r] = index[r];
      if (abs(sk[r]) == 1) //eqn (32) with algo4
      {
         beta = 1.0-alpha[r];
         ehat = (beta+beta*beta)*(0.5+alpha[r])*ehere+
                (alpha[r]+alpha[r]*alpha[r])*(0.5+beta)*ethere;//epsilon_hat
         rhere = ehere/ehat;//rho_here
         rthere = ethere/ehat;//rho_there
         a[0] = (beta+beta*beta)*rhere+alpha[r]*(1.0+2.0*beta)*rthere;//a_{-s}
         a[1] = -(beta+beta*beta)*rhere-(1.0+alpha[r])*(1.0+2.0*beta)*rthere;//a_{0}
         a[2] = (1.0+beta)*(1.0+beta)*rthere;//a_{s}
         a[3] = -beta*beta*rthere;//a_{2s}
         bk = -(beta+beta*beta)*(rthere-rhere);//b_k eqn(35)
         C = bk*tangent[r];
         for (s = 0; s < 4; s++)
         {
            tindex[r] = 2+(s-1)*sk[r];
            setvalarray(f[r],tindex,evalarray(f[r],tindex)+a[s]);
            //if sk = 1, f[r][1,2,3,4]coefficients of u[-1,0,1,2] from operator L eqn(32)
            //if sk = -1, f[r][0,1,2,3]coefficients of u[-2,-1,0,1] from operator L eqn(32)
         }
         tindex[r] = 2-sk[r];
         setvalarray(f[r],tindex,evalarray(f[r],tindex)-sk[r]*C*tangent[r]*sk[r]); //f[r][2-s] coeff of u[s] from operator T eqn (39) and eqn(40)
         tindex[r] = 2;
         setvalarray(f[r],tindex,evalarray(f[r],tindex)+sk[r]*C*tangent[r]*sk[r]);//f[r][2] coeff of u[0] from operator T eqn (39) and eqn(40)

         M[r][r] = 1.0-abs(sk[r])*(0.5+alpha[r])*bk*tangent[r]*tangent[r]; //eqn(40)

         for (j = 0; j < grid.dim; j++)
            if (j != r)
            {
               tmps = sk[j]; //tmps = s_j
               if (tmps == 0)
                  for (s = -1; s <= 1 && tmps == 0; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)//unused, tmp==0 already means no interface in j dim
                        tmps = s;
                  }
               if (tmps == 0)
               {
                  rindex[r] = index[r]-sk[r];
                  for (s = -1; s <= 1 && tmps == 0; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)//unused, tmp==0 already means no interface in j dim
                        tmps = s;
                  }
               }
               rindex[r] = index[r];
               rindex[j] = index[j];

               M[r][j] = -0.5*tmps*sk[r]*bk*tangent[j]*tangent[r]; //algo6, eqn(39)

               if (abs(tmps) == 1)
               {
                  tindex[j] = 2-tmps;
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                          sk[r]*C*tangent[j]*(1.0+alpha[r])*tmps); //coeff of u(x - sj h ej) in sk bk Tk
                  tindex[r] = 2-sk[r];
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                          sk[r]*C*tangent[j]*alpha[r]*tmps); //coeff of u(x - sk h ek - si h ej)
                  tindex[j] = 2;
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                          sk[r]*C*tangent[j]*alpha[r]*tmps); // coeff of u(x - sk h ek)
                  tindex[r] = 2;
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                          sk[r]*C*tangent[j]*(1.0+alpha[r])*tmps);//coeff of u(x)
               }
               else
               {
                  for (s = -1; s <= 1; s += 2)
                  {
                     tindex[j] = 2+s;
                     setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                             sk[r]*C*tangent[j]*(1.0+alpha[r])*0.5*s);
                  }
                  tindex[r] = 2-sk[r];
                  for (s = -1; s <= 1; s += 2)
                  {
                     tindex[j] = 2+s;
                     setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                             sk[r]*C*tangent[j]*alpha[r]*0.5*s);
                  }
                  tindex[r] = 2;
                  tindex[j] = 2;
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
            tindex[r] = 2+s;
            setvalarray(f[r],tindex,evalarray(f[r],tindex)+1.0);
         }
         tindex[r] = 2;
         setvalarray(f[r],tindex,evalarray(f[r],tindex)-2.0);
      }
   }
   // we get M uxx = D u + J, where M is 3x3, uxx is 3x1, uxx D is 3x(5x5x5) matrix, u is (5x5x5) vector
   gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);

   double uxx[grid.dim];
   forwardbacksub0(d,d,LU,PLR,PLC,grid.dim-1); //calculate inv(M)J
   for (r = 0; r < grid.dim; r++)
      uxx[r] = d[r]/(grid.dx[r]*grid.dx[r]);
    
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] < 5) //index is local index for u
   {
      for (r = 0; r < grid.dim; r++)
         temp[r] = evalarray(f[r],tindex);//temp is the 3x1 vector, one column of D, corresponding to u_tindex
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r]-2+tindex[r];//rindex is global index for u
      for (r = 0; r < grid.dim; r++)
         uxx[r] += temp[r]/(grid.dx[r]*grid.dx[r])*evalarray(u,rindex);//calculate inv(M)

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] >= 5; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
 
//   int m, sk2[grid.dim][grid.dim][4];
   int m, sk2[grid.dim][3][4];
   double tmperr[grid.dim]; 
//   double ***D1[grid.dim], ***D2[grid.dim][grid.dim], ***D1ucoef[grid.dim], **D1uxxcoef;
   double ***D1[grid.dim], ***D2[grid.dim][3], ***D1ucoef[grid.dim], **D1uxxcoef;
   for (r = 0; r < grid.dim; r++)
   {
      D1[r] = matrix(2,2,2);
      for (s = 0; s < grid.dim; s++)
         D2[r][s] = matrix(2,2,2);
   }
   for (r = 0; r < grid.dim; r++)
      D1ucoef[r] = matrix(2,2,2);
   D1uxxcoef = matrix(grid.dim-1,grid.dim-1);

   getD1(D1,sk,grid);
   for (r = 0; r < grid.dim; r++)
      for (s = r+1; s < grid.dim; s++)
         getsk2(sk2[r][s],r,s,index,S,grid);
   getD2(D2,sk2,grid);

   for (r = 0; r < grid.dim; r++)
      if (abs(sk[r]) == 1)
      {
         getDu(D1ucoef,D1uxxcoef,index,r,sk[r],alpha[r],thesign,sk,D1,D2,grid);
         for (s = 0; s < grid.dim; s++)
         {
            tmperr[s] = 0.0;
            for (m = 0; m < grid.dim; m++)
               tmperr[s] += D1uxxcoef[s][m]*uxx[m];
         }

         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] < 3)// go through neigbhor, calculate Du from D1ucoefficent * u
         {
            for (s = 0; s < grid.dim; s++)
               rindex[s] = index[s]-1+tindex[s];
            for (s = 0; s < grid.dim; s++)
               tmperr[s] += evalarray(D1ucoef[s],tindex)*evalarray(u,rindex);

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] >= 3; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }

         for (s = 0; s < grid.dim; s++)
         {
            tmperr[s] -= getDu(index,s,r,sk[r],alpha[r],thesign,grid);//getDu here calculate exact Du
            if (fabs(tmperr[s]) > fabs(uxerr[s]))
            {
               uxerr[s] = tmperr[s];
               globerrvec[s] = tmperr[s];
            }
            tmperr[s] = uxx[s]-getD2u(index,s,s,r,sk[r],alpha[r],thesign,grid);
            if (fabs(tmperr[s]) > fabs(globerrvec2[s]))
            {
               globerrvec2[s] = tmperr[s];
               for (m = 0; m < grid.dim; m++)
                  globerrmx[s][m] = index[m];
            }
         }
      }

   for (r = 0; r < grid.dim; r++)
   {
      free_matrix(D1[r],2,2,2);
      for (s = 0; s < grid.dim; s++)
         free_matrix(D2[r][s],2,2,2);
   }
   for (r = 0; r < grid.dim; r++)
      free_matrix(D1ucoef[r],2,2,2);
   free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);
   free_matrix(f,grid.dim-1,4,4,4);
}

void flipvector(double *v, double *w, int *index, int rstar, int sstar, double alpha,
                double *normal, double ***S, PBData &pb, GridData &grid)
{
// v can overwrite w
   int r, s, t;
   double ehere, ethere, thesign, sigma, Dudotn, Dtaudott, Dtaudots, length, 
          tol = 1.0e-14;
   double tangent1[grid.dim], tangent2[grid.dim], Dtau[grid.dim];

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

   length = 0.0;
   for (r = 0; r < grid.dim; r++)
   {   
      for (s = 0; s < grid.dim; s++)
         tangent2[s] = 0.0;
      tangent2[r] = 1.0;
      project(tangent2,normal,tangent2,grid.dim);
      if (sqrt(getdotprod(tangent2,tangent2,grid.dim)) > length)
      {
         for (t = 0; t < grid.dim; t++)
            tangent1[t] = tangent2[t];
         length = sqrt(getdotprod(tangent1,tangent1,grid.dim));
      }
   }
   for (s = 0; s < grid.dim; s++)
      tangent1[s] /= length;
   getunitcrossprod(tangent2,normal,tangent1);

   getDtau(Dtau,index,rstar,sstar,alpha,grid);
   getsigma(sigma,index,rstar,sstar,alpha,normal,pb,grid);
   Dudotn = getdotprod(w,normal,grid.dim);
   Dtaudott = getdotprod(Dtau,tangent1,grid.dim);
   Dtaudots = getdotprod(Dtau,tangent2,grid.dim);
   for (s = 0; s < grid.dim; s++)
      v[s] = w[s]-thesign*((sigma-thesign*(ehere-ethere)*Dudotn)*normal[s]/ethere+
                           Dtaudott*tangent1[s]+Dtaudots*tangent2[s]);
}

void flipvector2(double *v, double *w, int *index, int rstar, int sstar, double alpha,
                 double *normal, double ***S, PBData &pb, GridData &grid)
{
// v can overwrite w
   int r, s, t;
   double ehere, ethere, thesign, sigma, Dudotn, Dtaudott, Dtaudots, length, 
          tol = 1.0e-14;
   double tangent1[grid.dim], tangent2[grid.dim], Dtau[grid.dim], temp[grid.dim];

//   cout << "entering w " << w[0] << " " << w[1] << " " << w[2] << endl;

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

   length = 0.0;
   for (r = 0; r < grid.dim; r++)
   {   
      for (s = 0; s < grid.dim; s++)
         tangent2[s] = 0.0;
      tangent2[r] = 1.0;
      project(tangent2,normal,tangent2,grid.dim);
      if (sqrt(getdotprod(tangent2,tangent2,grid.dim)) > length)
      {
         for (t = 0; t < grid.dim; t++)
            tangent1[t] = tangent2[t];
         length = sqrt(getdotprod(tangent1,tangent1,grid.dim));
      }
   }
   for (s = 0; s < grid.dim; s++)
      tangent1[s] /= length;
   getunitcrossprod(tangent2,normal,tangent1);

   getDtau(Dtau,index,rstar,sstar,alpha,grid);
   getsigma(sigma,index,rstar,sstar,alpha,normal,pb,grid);
   Dudotn = getdotprod(w,normal,grid.dim);
   Dtaudott = getdotprod(Dtau,tangent1,grid.dim);
   Dtaudots = getdotprod(Dtau,tangent2,grid.dim);

/*
   double Duh[grid.dim], Dut[grid.dim], rhs[grid.dim], lhs[grid.dim];
   for (s = 0; s < grid.dim; s++)
   {
      Duh[s] = getDu(index,s,rstar,sstar,alpha,thesign,grid);
      Dut[s] = getDu(index,s,rstar,sstar,alpha,-thesign,grid);
   }
   temp[0] = getdotprod(normal,Dut,grid.dim);
   temp[1] = getdotprod(tangent1,Dut,grid.dim);
   temp[2] = getdotprod(tangent2,Dut,grid.dim);
   temp[0] *= ethere/ehere;
   for (s = 0; s < grid.dim; s++)
      lhs[s] = normal[s]*temp[0]+tangent1[s]*temp[1]+tangent2[s]*temp[2];
   cout << "real w " << Duh[0] << " " << Duh[1] << " " << Duh[2] << endl;
   cout << "approx w " << w[0] << " " << w[1] << " " << w[2] << endl;
   cout << "error " << fabs(w[0]-Duh[0]) << " " << fabs(w[1]-Duh[1]) 
        << " " << fabs(w[2]-Duh[2]) << endl;
   cout << "lhs " << lhs[0] << " " << lhs[1] << " " << lhs[2] << endl;
   cout << "lhs " << lhs[0] << " " << lhs[1] << " " << lhs[2] << endl;
   for (s = 0; s < grid.dim; s++)
      v[s] = Duh[s]-thesign*(sigma/ehere*normal[s]+Dtaudott*tangent1[s]+
                             Dtaudots*tangent2[s]);
   cout << "real rhs " << v[0] << " " << v[1] << " " << v[2] << endl;
*/
   for (s = 0; s < grid.dim; s++)
      v[s] = w[s]-thesign*(sigma/ehere*normal[s]+Dtaudott*tangent1[s]+
                           Dtaudots*tangent2[s]);
//   cout << "approx rhs " << v[0] << " " << v[1] << " " << v[2] << endl;

   temp[0] = getdotprod(normal,v,grid.dim);
   temp[1] = getdotprod(tangent1,v,grid.dim);
   temp[2] = getdotprod(tangent2,v,grid.dim);
   temp[0] /= ethere/ehere;
   for (s = 0; s < grid.dim; s++)
      v[s] = normal[s]*temp[0]+tangent1[s]*temp[1]+tangent2[s]*temp[2];
}

void getcim3Du(double &uint, double *Du, int *zindex, int rstar, int sstar, 
               double ***u, double ***S, PBData &pb, GridData &grid)
{
// switches to computations on other side if current point not cim3
// stencil = 2
   int r, s, t, i, m, n;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim], index[grid.dim];
   double d0[grid.dim];
   double ***dcoef[grid.dim];
//   double ***D1[grid.dim], ***D2[grid.dim][grid.dim];
   double ***D1[grid.dim], ***D2[grid.dim][3];
   double **LU, **G;
//   int gamma[grid.dim][2], sk2[grid.dim][grid.dim][4];
   int gamma[grid.dim][2], sk2[grid.dim][3][4];
   int PLR[grid.dim], PLC[grid.dim];
   double alpha, beta, finalalpha, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[grid.dim];
   int sk[grid.dim];
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;
   double jumpD1u, ***jumpD1ucoef, jumpD1uxxcoef[grid.dim];
   double jumpD2u, ***jumpD2ucoef, jumpD2uxxcoef[grid.dim];
   double ****D1ucoef[grid.dim], ***D1uxxcoef;
   double tempr1, tempr2, sigmar; 
   double ux, uxx[grid.dim];

   for (r = 0; r < grid.dim; r++)
      index[r] = zindex[r];
   if (checkcim3(S,index,grid) == 3)
   {
      getcim4Du(uint,Du,index,rstar,sstar,u,S,pb,grid);
      cout << "cim 4 Du = " << Du[0] << " " << Du[1] << " " << Du[2] << endl;
      index[rstar] = zindex[rstar]+sstar;
      if (evalarray(S,index) < 0.0)
         thesign = -1.0;
      else
         thesign = 1.0;
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
      for (r = 0; r < grid.dim; r++)
         sk[r] = gamma[r][1]-gamma[r][0];
      getinterfaceinfo(alpha,tangent,normal,S,index,zindex,grid);
      cout << "exact " << getDu(index,0,rstar,sk[rstar],alpha,-thesign,grid) << " "
           << getDu(index,1,rstar,sk[rstar],alpha,-thesign,grid) << " "
           << getDu(index,2,rstar,sk[rstar],alpha,-thesign,grid) << endl;
      if (checkcim3(S,index,grid) == 2)
      {
         getcim3Du(uint,Du,index,rstar,sk[rstar],u,S,pb,grid);
         flipvector(Du,Du,index,rstar,sk[rstar],alpha,normal,S,pb,grid);
         gettau(tau,index,r,sk[rstar],alpha,grid);
         uint -= thesign*tau;
         cout << "flip " << Du[0] << " " << Du[1] << " " << Du[2] << endl;
//         double realDu[grid.dim], junk[grid.dim];
//         cout << "index " << index[0] << " " << index[1] << " " << index[2] << endl;
//         cout << "status " << checkcim3(S,index,grid) << endl;
//         cout << "calc Du " << Du[0] << " " << Du[1] << " " << Du[2] << endl;
//         flipvector(junk,Du,index,rstar,sk[rstar],alpha,normal,S,pb,grid);
//         cout << "approx flip 1 " << junk[0] << " " << junk[1] << " " << junk[2] << endl;
//         flipvector2(junk,Du,index,rstar,sk[rstar],alpha,normal,S,pb,grid);
//         cout << "approx flip 2 " << junk[0] << " " << junk[1] << " " << junk[2] << endl;
//         for (s = 0; s < grid.dim; s++)
//            realDu[s] = getDu(index,s,rstar,sk[rstar],alpha,thesign,grid);
//         flipvector(junk,realDu,index,rstar,sk[rstar],alpha,normal,S,pb,grid);
//         cout << "exact flip 1 " << junk[0] << " " << junk[1] << " " << junk[2] << endl;
//         flipvector2(junk,realDu,index,rstar,sk[rstar],alpha,normal,S,pb,grid);
//         cout << "exact flip 2 " << junk[0] << " " << junk[1] << " " << junk[2] << endl;
      }
      else
      {
         cout << "some kind of error" << endl;
         cout << index[0] << " " << index[1] << " " << index[2] << endl;
         cout << checkcim3(S,rindex,grid) << endl;
         for (r = 0; r < grid.dim; r++)
            rindex[r] = index[r];
         for (r = 0; r < grid.dim; r++)
         {
            for (s = -1; s <= 1; s += 2)
            {
               rindex[r] = index[r]+s;
               cout << checkcim3(S,rindex,grid) << " ";
            }
            rindex[r] = index[r];
            cout << endl;
         }
         cout << rstar << " " << sstar << endl;
         exit(1);
      }
      getchar();
   }
   else
   {
      LU = matrix(grid.dim-1,grid.dim-1);
      G = matrix(grid.dim-1,grid.dim-1);
      jumpD1ucoef = matrix(2,2,2);
      jumpD2ucoef = matrix(2,2,2);
      for (r = 0; r < grid.dim; r++)
      {
         D1[r] = matrix(2,2,2);
         for (s = 0; s < grid.dim; s++)
            D2[r][s] = matrix(2,2,2);
         dcoef[r] = matrix(2,2,2);
      }
      for (r = 0; r < grid.dim; r++)
         D1ucoef[r] = matrix(grid.dim-1,2,2,2);
      D1uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   
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
      for (r = 0; r < grid.dim; r++)
         sk[r] = gamma[r][1]-gamma[r][0];
// get sk2
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
            getsk2(sk2[m][n],m,n,index,S,grid);
   
      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 1;
      for (r = 0; r < grid.dim; r++)
         if (abs(sk[r]) == 1)
         {
            rindex[r] = index[r]+sk[r];
// getting first and second derivatives
            getD1(D1,sk,grid);
            getD2(D2,sk2,grid);
// getting interface info
            getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
            if (r == rstar && sk[r] == sstar)
               finalalpha = alpha;
            beta = 1.0-alpha;
// get Du 2nd order
            getDu(D1ucoef[r],D1uxxcoef[r],index,r,sk[r],alpha,thesign,sk,D1,D2,grid);
// getting jump of uxx in r
            gettau(tau,index,r,sk[r],alpha,grid);
            getjumpux(jumpD1u,jumpD1ucoef,jumpD1uxxcoef,index,r,sk,alpha,
                      thesign,normal,tangent,D1ucoef[r],D1uxxcoef[r],S,pb,grid);
            if (globorder == 1)
               getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,index,r,sk,alpha,thesign,
                          normal,D1,D2,S,pb,grid);
            else
               getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxxcoef,index,r,sk,alpha,thesign,
                          normal,D1ucoef[r],D1uxxcoef[r],D2,S,pb,grid);
// form d0 and dcoef's rth entry 
            d0[r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+sk[r]*beta*jumpD1u/grid.dx[r]+
                             0.5*beta*beta*jumpD2u);
   
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] < 3)
            {
               setvalarray(dcoef[r],sindex,
                           thesign*(sk[r]*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                    0.5*beta*beta*evalarray(jumpD2ucoef,sindex))-
                           sk[r]*(alpha+beta)*evalarray(D1ucoef[r][r],sindex)/
                           grid.dx[r]);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 1;
            setvalarray(dcoef[r],sindex,
                        evalarray(dcoef[r],sindex)-1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 1+sk[r];
            setvalarray(dcoef[r],sindex,
                        evalarray(dcoef[r],sindex)+1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 1;
// form matrix's rth row
            for (m = 0; m < grid.dim; m++)
               G[r][m] = -thesign*(sk[r]*beta*jumpD1uxxcoef[m]/grid.dx[r]+
                                   0.5*beta*beta*jumpD2uxxcoef[m])+
                         sk[r]*(alpha+beta)*D1uxxcoef[r][r][m]/grid.dx[r];
            G[r][r] += 0.5*(beta*beta-alpha*alpha);
   
            rindex[r] = index[r];
         }
         else
         {
// if interior point
            d0[r] = 0.0;
            for (s = 0; s < grid.dim; s++)
               G[r][s] = 0.0;
            G[r][r] = 1.0;
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] < 3)
            {
               setvalarray(dcoef[r],sindex,0.0);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 1;
            setvalarray(dcoef[r],sindex,-2.0/(grid.dx[r]*grid.dx[r]));
            for (s = -1; s <= 1; s += 2)
            {
               sindex[r] = 1+s;
               setvalarray(dcoef[r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
            }
            sindex[r] = 1;
         }
   
      gecp0(LU,PLR,PLC,G,grid.dim-1,grid.dim-1);
      forwardbacksub0(temp,d0,LU,PLR,PLC,grid.dim-1);
      for (s = 0; s < grid.dim; s++)
         uxx[s] = temp[s];
      ux = 0.0;
   
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s];
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 0;
      while (sindex[0] < 3)
      {
         for (n = 0; n < grid.dim; n++)
            temp[n] = evalarray(dcoef[n],sindex);
         forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
         for (s = 0; s < grid.dim; s++)
            tindex[s] = index[s]-1+sindex[s];
         for (s = 0; s < grid.dim; s++)
            uxx[s] += temp[s]*evalarray(u,tindex);
         ux += evalarray(D1[rstar],sindex)*evalarray(u,tindex);

         (sindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
         {
            sindex[i] = 0;
            (sindex[i-1])++;
         }
      }
      for (s = 0; s < grid.dim; s++)
         sindex[s] = 1;
   
      r = rstar;
      if (sk[r] == sstar)
      {
         if (globintorder == 3)
         {
            ux += sk[r]*0.5*grid.dx[r]*uxx[r];
            uint = evalarray(u,index)+sk[r]*finalalpha*grid.dx[r]*ux+
                   0.5*(finalalpha*grid.dx[r])*(finalalpha*grid.dx[r])*uxx[r];
         }
         else
            uint = evalarray(u,index)+sk[r]*finalalpha*grid.dx[r]*ux;

         for (s = 0; s < grid.dim; s++)
         {
            Du[s] = 0.0;
            for (m = 0; m < grid.dim; m++)
               Du[s] += D1uxxcoef[r][s][m]*uxx[m];
         }

         for (s = 0; s < grid.dim; s++)
            sindex[s] = 0;
         while (sindex[0] < 3)
         {
            for (s = 0; s < grid.dim; s++)
               tindex[s] = index[s]-1+sindex[s];
            for (s = 0; s < grid.dim; s++)
               Du[s] += evalarray(D1ucoef[r][s],sindex)*evalarray(u,tindex);
   
            (sindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && sindex[i] >= 3; i--)
            {
               sindex[i] = 0;
               (sindex[i-1])++;
            }
         }
         for (s = 0; s < grid.dim; s++)
            sindex[s] = 1;
      }

      free_matrix(LU,grid.dim-1,grid.dim-1);
      free_matrix(G,grid.dim-1,grid.dim-1);
      free_matrix(jumpD1ucoef,2,2,2);
      free_matrix(jumpD2ucoef,2,2,2);
      for (r = 0; r < grid.dim; r++)
      {
         free_matrix(D1[r],2,2,2);
         for (s = 0; s < grid.dim; s++)
            free_matrix(D2[r][s],2,2,2);
         free_matrix(dcoef[r],2,2,2);
      }
      for (r = 0; r < grid.dim; r++)
         free_matrix(D1ucoef[r],grid.dim-1,2,2,2);
      free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   }
}

void checkcim3Du(double ***u, double ***S, PBData &pb, GridData &grid)
{
   int i, r, s, t, m, thestatus; 
   int tindex[grid.dim], rindex[grid.dim], zindex[grid.dim][grid.dim];
   double temp, uint, Du[grid.dim], theerr[grid.dim];
   double alpha, thesign, tangent[grid.dim], normal[grid.dim];
   int eindex[grid.dim];
   eindex[0] = 25;
   eindex[1] = 34;
   eindex[2] = 40;

   for (i = 0; i < grid.dim; i++)
      theerr[i] = 0.0;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      thestatus = checkcim3(S,tindex,grid);
/*
      if ((tindex[0] == 25 && tindex[1] == 34 && tindex[2] == 40) ||
          (tindex[0] == 25 && tindex[1] == 39 && tindex[2] == 17) ||
          (tindex[0] == 25 && tindex[1] == 21 && tindex[2] == 7))
      {
         cout << "in checkcim3Du, status = " << thestatus << endl;
         for (r = 0; r < grid.dim; r++)
            rindex[r] = tindex[r];
         for (r = 0; r < grid.dim; r++)
         {
            for (s = -1; s <= 1; s += 2)
            {
               rindex[r] = tindex[r]+s;
               cout << "   neighboring status = " << checkcim3(S,rindex,grid) << endl;
            }
            rindex[r] = tindex[r];
         }
      }
*/
      if (thestatus == 2 || thestatus == 3)
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
                  getcim3Du(uint,Du,tindex,r,s,u,S,pb,grid);
                  cout << "point status = " << checkcim3(S,tindex,grid) << endl;
                  cout << "unflipped" << endl;
                  cout << Du[0] << " " << Du[1] << " " << Du[2] << endl;
                  getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
                  cout << getDu(tindex,0,r,s,alpha,thesign,grid) << " " 
                       << getDu(tindex,1,r,s,alpha,thesign,grid) << " " 
                       << getDu(tindex,2,r,s,alpha,thesign,grid) << endl;
                  cout << fabs(Du[0]-getDu(tindex,0,r,s,alpha,thesign,grid)) << " "
                       << fabs(Du[1]-getDu(tindex,1,r,s,alpha,thesign,grid)) << " "
                       << fabs(Du[2]-getDu(tindex,2,r,s,alpha,thesign,grid)) << endl;
                  int aindex[grid.dim];
                  double junk[grid.dim];
                  for (t = 0; t < grid.dim; t++)
                     aindex[t] = tindex[t];
                  aindex[r] = tindex[r]+s;
                  flipvector(junk,Du,tindex,r,s,alpha,normal,S,pb,grid);
                  cout << "point status = " << checkcim3(S,aindex,grid) << endl;
                  cout << "flipped" << endl;
                  cout << junk[0] << " " << junk[1] << " " << junk[2] << endl;
                  cout << getDu(tindex,0,r,s,alpha,-thesign,grid) << " " 
                       << getDu(tindex,1,r,s,alpha,-thesign,grid) << " "
                       << getDu(tindex,2,r,s,alpha,-thesign,grid) << endl;
                  cout << fabs(junk[0]-getDu(tindex,0,r,s,alpha,-thesign,grid)) << " "
                       << fabs(junk[1]-getDu(tindex,1,r,s,alpha,-thesign,grid)) << " "
                       << fabs(junk[2]-getDu(tindex,2,r,s,alpha,-thesign,grid)) << endl;
                  cout << "flipped 2" << endl;
                  flipvector2(junk,Du,tindex,r,s,alpha,normal,S,pb,grid);
                  cout << junk[0] << " " << junk[1] << " " << junk[2] << endl;
                  cout << getDu(tindex,0,r,s,alpha,-thesign,grid) << " " 
                       << getDu(tindex,1,r,s,alpha,-thesign,grid) << " "
                       << getDu(tindex,2,r,s,alpha,-thesign,grid) << endl;
                  cout << fabs(junk[0]-getDu(tindex,0,r,s,alpha,-thesign,grid)) << " "
                       << fabs(junk[1]-getDu(tindex,1,r,s,alpha,-thesign,grid)) << " "
                       << fabs(junk[2]-getDu(tindex,2,r,s,alpha,-thesign,grid)) << endl;
                  cout << "computed flipped" << endl;
                  getcim3Du(uint,junk,aindex,r,-s,u,S,pb,grid);
                  cout << junk[0] << " " << junk[1] << " " << junk[2] << endl;
                  cout << getDu(tindex,0,r,s,alpha,-thesign,grid) << " " 
                       << getDu(tindex,1,r,s,alpha,-thesign,grid) << " "
                       << getDu(tindex,2,r,s,alpha,-thesign,grid) << endl;
                  cout << fabs(junk[0]-getDu(tindex,0,r,s,alpha,-thesign,grid)) << " "
                       << fabs(junk[1]-getDu(tindex,1,r,s,alpha,-thesign,grid)) << " "
                       << fabs(junk[2]-getDu(tindex,2,r,s,alpha,-thesign,grid)) << endl;
                  getchar();
                  for (t = 0; t < grid.dim; t++)
                  {
                     temp = getDu(tindex,t,r,s,alpha,thesign,grid);
                     if (fabs(theerr[t]) < fabs(temp-Du[t]))
                     {
                        theerr[t] = temp-Du[t];
                        for (m = 0; m < grid.dim; m++)  
                           zindex[t][m] = tindex[m];
                     }
                  }

               }
            }
            rindex[r] = tindex[r];
         }
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "cim3 Du error = " << theerr[0] << " " << theerr[1] << " " << theerr[2] 
        << endl;
   for (s = 0; s < grid.dim; s++)
      cout << "at " << zindex[s][0] << " " << zindex[s][1] << " " << zindex[s][2] 
           << endl;
}

void checkallDu(double ***u, double ***S, PBData &pb, GridData &grid)
{
   int i, r, s, t;
   int tindex[grid.dim], rindex[grid.dim];
   double uint, Du[grid.dim];
   double alpha, thesign, tangent[grid.dim], normal[grid.dim];
   int thenum, two23, flip, thestatus[2*grid.dim+1], count[12], catcount[4], 
       betcount[6][2], anocount[2][12], anocat[2][4];
   double tmperr, theerr[12][grid.dim], maxerr[grid.dim], temp[grid.dim], 
          anoerr[2][12][grid.dim], anomax[2][grid.dim];
   for (r = 0; r < 12; r++)
   {
      count[r] = 0;
      anocount[0][r] = 0;
      anocount[1][r] = 0;
      for (s = 0; s < grid.dim; s++)
      {
         theerr[r][s] = 0.0;
         anoerr[0][r][s] = 0.0;
         anoerr[1][r][s] = 0.0;
      }
   }
   for (r = 0; r < grid.dim; r++)
   {
      maxerr[r] = 0.0;
      anomax[0][r] = 0.0;
      anomax[1][r] = 0.0;
   }

   for (r = 0; r < 4; r++)
   {
      catcount[r] = 0;
      anocat[0][r] = 0;
      anocat[1][r] = 0;
   }
   for (r = 0; r < 6; r++)
      for (s = 0; s < 2; s++)
         betcount[r][s] = 0;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (globcim == 4)
         thestatus[0] = checkcim3(S,tindex,grid);
      else if (globcim == 5 || globcim == 345)
         thestatus[0] = getstatus5(S,tindex,grid);
      if (thestatus[0] == 4)
      {
         cout << "Found a status 4 point" << endl;
         exit(1);
      }
      else if (thestatus[0] == 2 || thestatus[0] == 3)
      {
         for (r = 0; r < grid.dim; r++)
            rindex[r] = tindex[r];
         for (r = 0; r < grid.dim; r++)
         {
            for (s = -1; s <= 1; s += 2)
            {
               rindex[r] = tindex[r]+s;
               if (globcim == 4)
                  thestatus[2*r+(s+3)/2] = checkcim3(S,rindex,grid);
               else if (globcim == 5 || globcim == 345)
                  thestatus[2*r+(s+3)/2] = getstatus5(S,rindex,grid);
            }
            rindex[r] = tindex[r];
         }
         for (r = 1; r < 2*grid.dim+1 && thestatus[r] <= 2; r++);
         if (thestatus[0] == 2 && r >= 2*grid.dim+1)
         {
            thenum = 0;
            (catcount[0])++;
            if (thesign < 0.0)
               (anocat[0][0])++;
            else
               (anocat[1][0])++;
         }
         else if (thestatus[0] == 2)
         {
            thenum = 2;
            (catcount[1])++;
            if (thesign < 0.0)
               (anocat[0][1])++;
            else
               (anocat[1][1])++;
         }
         else if (thestatus[0] == 3 && r >= 2*grid.dim+1)
         {
            thenum = 6;
            (catcount[2])++;
            if (thesign < 0.0)
               (anocat[0][2])++;
            else
               (anocat[1][2])++;
         }
         else if (thestatus[0] == 3)
         {
            thenum = 8;
            (catcount[3])++;
            if (thesign < 0.0)
               (anocat[0][3])++;
            else
               (anocat[1][3])++;
            cout << tindex[0] << " " << tindex[1] << " " << tindex[2] << " " 
                 << evalarray(S,tindex) << endl;
         }
         else
         {
            cout << "Missed some case ";
            for (r = 0; r < 2*grid.dim+1; r++)
               cout << thestatus[r] << " ";
            cout << endl;
            exit(1);
         }

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
               if ((evalarray(S,tindex) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)
               {
                  if (thenum == 0 || thenum == 2)
                  {
                     if (!globallcim345)
                        getcim3Du(uint,Du,tindex,r,s,u,S,pb,grid);
                     else
                        getcim345Du(uint,Du,tindex,r,s,u,S,pb,grid);
                  }
                  else if (thenum == 6 || thenum == 8)
                  {
                     if (globcim == 4)
                        getcim4Du(uint,Du,tindex,r,s,u,S,pb,grid);
                     else if (globcim == 5)
                        getcim5Du(uint,Du,tindex,r,s,u,S,pb,grid);
                     else if (globcim == 345)
                        getcim345Du(uint,Du,tindex,r,s,u,S,pb,grid);
                  }
                  getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
                  if (thestatus[2*r+(s+3)/2] == 2)
                     two23 = 0;
                  else
                     two23 = 1;
                  flip = 0;
                  (count[thenum+flip+2*two23])++;
                  if (thesign < 0.0)
                     (anocount[0][thenum+flip+2*two23])++;
                  else
                     (anocount[1][thenum+flip+2*two23])++;
                  for (t = 0; t < grid.dim; t++)
                  {
                     tmperr = fabs(Du[t]-getDu(tindex,t,r,s,alpha,thesign,grid));
                     if (tmperr > theerr[thenum+flip+2*two23][t])
                        theerr[thenum+flip+2*two23][t] = tmperr;
                     if (thesign < 0.0 && tmperr > anoerr[0][thenum+flip+2*two23][t])
                        anoerr[0][thenum+flip+2*two23][t] = tmperr;
                     else if (thesign >= 0.0 && 
                              tmperr > anoerr[1][thenum+flip+2*two23][t])
                        anoerr[1][thenum+flip+2*two23][t] = tmperr;
                  }
                  flipvector(Du,Du,tindex,r,s,alpha,normal,S,pb,grid);
                  flip = 1;
                  (count[thenum+flip+2*two23])++;
                  if (thesign < 0.0)
                     (anocount[0][thenum+flip+2*two23])++;
                  else
                     (anocount[1][thenum+flip+2*two23])++;
                  for (t = 0; t < grid.dim; t++)
                  {
                     tmperr = fabs(Du[t]-getDu(tindex,t,r,s,alpha,-thesign,grid));
                     if (tmperr > theerr[thenum+flip+2*two23][t])
                        theerr[thenum+flip+2*two23][t] = tmperr;
                     if (thesign < 0.0 && tmperr > anoerr[0][thenum+flip+2*two23][t])
                        anoerr[0][thenum+flip+2*two23][t] = tmperr;
                     else if (thesign >= 0.0 && 
                              tmperr > anoerr[1][thenum+flip+2*two23][t])
                        anoerr[1][thenum+flip+2*two23][t] = tmperr;
                  }
/*
                  if (thenum == 0 || thenum == 2)
                  {
                     getcim3Du(temp,rindex,r,-s,u,S,pb,grid);
                     for (t = 0; t < grid.dim; t++)
                        if (fabs(Du[t]-getDu(tindex,t,r,s,alpha,-thesign,grid)) <
                            fabs(temp[t]-getDu(tindex,t,r,s,alpha,-thesign,grid)))
                           (betcount[thenum+flip+2*two23][0])++;
                  }
*/
               }
            }
            rindex[r] = tindex[r];
         }
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   for (r = 0; r < 12; r += 2)
      for (s = 0; s < grid.dim; s++)
      {
         if (theerr[r][s] > maxerr[s])
            maxerr[s] = theerr[r][s];
         if (anoerr[0][r][s] > anomax[0][s])
            anomax[0][s] = anoerr[0][r][s];
         if (anoerr[1][r][s] > anomax[1][s])
            anomax[1][s] = anoerr[1][r][s];
      }

   cout << "Total count: ";
   for (r = 0; r < 4; r++)
      cout << catcount[r] << " ";
   cout << endl;
   cout << "Count and Errors:" << endl;
   for (r = 0; r < 12; r++)
   {
      cout << "   " << count[r] << " ";
      for (s = 0; s < grid.dim; s++)
         cout << theerr[r][s] << " ";
      cout << endl;
   }
   cout << "Direct Error:" << endl;
   for (r = 0; r < grid.dim; r++)
      cout << maxerr[r] << " ";
   cout << endl;
   cout << "Downloadable:" << endl;
   cout << grid.nx[0] << " ";
   for (r = 0; r < 4; r++)
      cout << catcount[r] << " ";
   for (r = 0; r < 12; r++)
      for (s = 0; s < grid.dim; s++)
         cout << theerr[r][s] << " ";
   for (r = 0; r < 12; r++)
      cout << count[r] << " ";
   for (r = 0; r < grid.dim; r++)
      cout << maxerr[r] << " ";
   cout << endl;
   cout << "Downloadable -:" << endl;
   cout << grid.nx[0] << " ";
   for (r = 0; r < 4; r++)
      cout << anocat[0][r] << " ";
   for (r = 0; r < 12; r++)
      for (s = 0; s < grid.dim; s++)
         cout << anoerr[0][r][s] << " ";
   for (r = 0; r < 12; r++)
      cout << anocount[0][r] << " ";
   for (r = 0; r < grid.dim; r++)
      cout << anomax[0][r] << " ";
   cout << endl;
   cout << "Downloadable +:" << endl;
   cout << grid.nx[0] << " ";
   for (r = 0; r < 4; r++)
      cout << anocat[1][r] << " ";
   for (r = 0; r < 12; r++)
      for (s = 0; s < grid.dim; s++)
         cout << anoerr[1][r][s] << " ";
   for (r = 0; r < 12; r++)
      cout << anocount[1][r] << " ";
   for (r = 0; r < grid.dim; r++)
      cout << anomax[1][r] << " ";
   cout << endl;
}

void checkcim5Du(double ***u, double ***S, PBData &pb, GridData &grid)
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
      thestatus = getstatus5(S,tindex,grid);
      if (thestatus == 4)
         cout << "CIM 4 point at " << tindex[0] << " " << tindex[1] << " " 
              << tindex[2] << endl;
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
               if (thestatus == 2)
                  getcim3Du(uint,Du,tindex,r,s,u,S,pb,grid);
               else if (thestatus == 3)
                  getcim5Du(uint,Du,tindex,r,s,u,S,pb,grid);
               else if (thestatus == 4)
                  getcim4Du(uint,Du,tindex,r,s,u,S,pb,grid);
               else
                  cout << "bad for status " << thestatus << " at " 
                       << tindex[0] << " " << tindex[1] << " " << tindex[2] 
                       << " in direction " << r << " " << s << endl;

               getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
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

   cout << "Error in cim5 derivative is " << theerr[0] << " " << theerr[1] << " " 
        << theerr[2] << endl;
   cout << "   at";
   for (t = 0; t < grid.dim; t++)
   {
      for (m = 0; m < grid.dim; m++)
         cout << " " << zindex[t][m];
      cout << " (" << getstatus5(S,zindex[t],grid) << ")";
      if (t < grid.dim-1)
         cout << ", ";
      else
         cout << "." << endl;
   }

   free_matrix(zindex,grid.dim-1,grid.dim-1);
}
// checkdu NO StorageStruct
void checkcim345Du(double ***u, double ***S, PBData &pb, GridData &grid)
{
   cout<< "[checkcim345Du]" <<endl;
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
      thestatus = getstatus5(S,tindex,grid); 
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
               if (thestatus == 2){
                  if (!globallcim345)
                     getcim3Du(uint,Du,tindex,r,s,u,S,pb,grid);
                  else
                     getcim345Du(uint,Du,tindex,r,s,u,S,pb,grid);
               }else if (thestatus == 3 || thestatus == 4){
                  getcim345Du(uint,Du,tindex,r,s,u,S,pb,grid);
               }else{
                  getcim12Du(uint,Du,tindex,r,s,u,S,pb,grid);
                  print_surf(S,tindex,2);
                  cout << "bad for status " << thestatus << " at " << tindex[0] << " " << tindex[1] << " " << tindex[2] << endl;
               }
//               getinterfaceDu(uint,Du,tindex,r,s,u,S,pb,grid);
               getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
               for (t = 0; t < grid.dim; t++)
               {
                  tmperr = fabs(Du[t]-getDu(tindex,t,r,s,alpha,thesign,grid));//error in dim t
                  if (tmperr > theerr[t])
                  {
                     theerr[t] = tmperr; //theerr[t] is max error in dim t
                     for (m = 0; m < grid.dim; m++)
                        zindex[t][m] = tindex[m]; //zindex[t] is index of max err in dim t
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
   
   cout << "Error in cim345 derivative is " << theerr[0] << " " << theerr[1] << " " 
        << theerr[2] << endl;
   cout << "   at";
   for (t = 0; t < grid.dim; t++)
   {
      for (m = 0; m < grid.dim; m++)
         cout << " " << zindex[t][m];
      cout << " (" << getstatus5(S,zindex[t],grid) << ")";
      if (t < grid.dim-1)
         cout << ", ";
      else
         cout << "." << endl;
   }

   free_matrix(zindex,grid.dim-1,grid.dim-1);
}
// check du with StorageStruct at interface
void checkcim345Du(double ***u, StorageStruct *Dusmall, int smallsize, double ***S, 
                   PBData &pb, GridData &grid)
{
  cout<< "[checkcim345Du Storage]" <<endl;
   int i, r, s, t, m, mid = 2;
   int tindex[grid.dim], rindex[grid.dim], **zindex = imatrix(grid.dim-1,grid.dim-1);
   double uint, Du[grid.dim];
   double alpha, thesign, tangent[grid.dim], normal[grid.dim];
   double tmperr, theerr[grid.dim];
   for (t = 0; t < grid.dim; t++)
      theerr[t] = 0.0;

   double uint2, Du2[grid.dim];

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
//               getinterfaceDu(uint2,Du2,tindex,r,s,u,S,pb,grid);
//               evalfromstorage(uint,Du,tindex,r,s,mid,Dusmall,smallsize,u,S,pb,grid);
               evalfromstorage(uint,Du,tindex,r,s,Dusmall,smallsize,u,S,pb,grid);

               getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
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

               if (globwriteerr){
                  outfile_Duerr<<tindex[0]<<","<<tindex[1]<<","<<tindex[2]<<","<<r<<","<<s<<",";
                  outfile_Duerr <<setprecision(12)
                  <<Du[0]-getDu(tindex,0,r,s,alpha,thesign,grid)<<","
                  <<Du[1]-getDu(tindex,1,r,s,alpha,thesign,grid)<<","
                  <<Du[2]-getDu(tindex,2,r,s,alpha,thesign,grid)
                  <<endl;;
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
   
   cout << "Error in cim345 derivative is " << theerr[0] << " " << theerr[1] << " " 
        << theerr[2] << endl;
   cout << "   at";
   for (t = 0; t < grid.dim; t++)
   {
      for (m = 0; m < grid.dim; m++)
         cout << " " << zindex[t][m];
      cout << " (" << getstatus5(S,zindex[t],grid) << ")";
      if (t < grid.dim-1)
         cout << ", ";
      else
         cout << "." << endl;
   }

   free_matrix(zindex,grid.dim-1,grid.dim-1);
}
// check cim 345 Du in positive region and negative region
void checkcim345Dupm(double ***u, double ***S, PBData &pb, GridData &grid)
{
  cout << "[checkcim345Dupm]"<<endl;
   int i, r, s, t, m, pm, thestatus;
   int tindex[grid.dim], rindex[grid.dim], **zindex[2], **vindex, **windex;
   double uint[2], **Du = matrix(1,grid.dim-1);
   double alpha, thesign, tangent[grid.dim], normal[grid.dim];
   double tmperr, Dtheerr[2][grid.dim], theerr[2], valerr[2];

   for (s = 0; s < 2; s++)
   {
      theerr[s] = 0.0;
      for (t = 0; t < grid.dim; t++)
         Dtheerr[s][t] = 0.0;
      valerr[s] = 0.0;
   }
   vindex = imatrix(1,grid.dim-1);
   for (s = 0; s < 2; s++)
      zindex[s] = imatrix(grid.dim-1,grid.dim-1);
   windex = imatrix(1,grid.dim-1);

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      thestatus = getstatus5(S,tindex,grid);
      if (evalarray(S,tindex) < 0.0)
      {
         thesign = -1.0;
         pm = 0;
      }
      else
      {
         thesign = 1.0;
         pm = 1;
      }
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
               if (thestatus == 2)
                  if (!globallcim345)
                     getcim3Du(uint[pm],Du[pm],tindex,r,s,u,S,pb,grid);
                  else
                     getcim345Du(uint[pm],Du[pm],tindex,r,s,u,S,pb,grid);
               else if (thestatus == 3 || thestatus == 4)
                  getcim345Du(uint[pm],Du[pm],tindex,r,s,u,S,pb,grid);
               else{
                   getcim12Du(uint[pm],Du[pm],tindex,r,s,u,S,pb,grid);
                   print_surf(S,tindex,2);
                   cout << "bad for status " << thestatus << " at " << tindex[0] << " " << tindex[1] << " " << tindex[2] << endl;
                }

               tmperr = fabs(evalarray(u,tindex)-getu(tindex,0,0,0.0,thesign,grid));//u error at grid point
               if (tmperr > theerr[pm])
               {
                  theerr[pm] = tmperr;
                  for (m = 0; m < grid.dim; m++)
                     vindex[pm][m] = tindex[m];
               }

               getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
               for (t = 0; t < grid.dim; t++)
               {
                  tmperr = fabs(Du[pm][t]-getDu(tindex,t,r,s,alpha,thesign,grid));//Du error at interface
                  if (tmperr > Dtheerr[pm][t])
                  {
                     Dtheerr[pm][t] = tmperr;
                     for (m = 0; m < grid.dim; m++)
                        zindex[pm][t][m] = tindex[m];
                  }
               }

               tmperr = fabs(uint[pm]-getu(tindex,r,s,alpha,thesign,grid)); // u error at interface
               if (tmperr > valerr[pm])
               {
                  valerr[pm] = tmperr;
                  
                  for (m = 0; m < grid.dim; m++)
                     windex[pm][m] = tindex[m];
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
   

   cout << "Error in cim345 in negative region is " << theerr[0] << endl;
   cout << "   at " << vindex[0][0] << " " << vindex[0][1] << " " << vindex[0][2] 
        << " (" << getstatus5(S,vindex[0],grid) << ")." << endl<<endl;

   cout << "Error in cim345 interface in negative region is " << valerr[0] << endl;
   cout << "   at " << windex[0][0] << " " << windex[0][1] << " " << windex[0][2] 
        << " (" << getstatus5(S,windex[0],grid) << ")." << endl<<endl;

   cout << "Error in cim345 derivative in negative region is " 
        << Dtheerr[0][0] << " " << Dtheerr[0][1] << " " << Dtheerr[0][2] << endl;
        
   
   cout << "   at";
   for (t = 0; t < grid.dim; t++)
   {
      for (m = 0; m < grid.dim; m++)
         cout << " " << zindex[0][t][m];
      cout << " (" << getstatus5(S,zindex[0][t],grid) << ")";
      if (t < grid.dim-1)
         cout << ", ";
      else
         cout << "." << endl;
   }
   cout << "Error in cim345 in positive region is " << theerr[1] << endl;
   cout << "   at " << vindex[1][0] << " " << vindex[1][1] << " " << vindex[1][2] 
        << " (" << getstatus5(S,vindex[1],grid) << ")." << endl<<endl;

   cout << "Error in cim345 interface in positive region is " << valerr[1] << endl;
   cout << "   at " << windex[1][0] << " " << windex[1][1] << " " << windex[1][2] 
        << " (" << getstatus5(S,windex[1],grid) << ")." << endl<<endl;
   
   cout << "Error in cim345 derivative in positive region is " 
        << Dtheerr[1][0] << " " << Dtheerr[1][1] << " " << Dtheerr[1][2] << endl;
   cout << "   at";
   for (t = 0; t < grid.dim; t++)
   {
      for (m = 0; m < grid.dim; m++)
         cout << " " << zindex[1][t][m];
      cout << " (" << getstatus5(S,zindex[1][t],grid) << ")";
      if (t < grid.dim-1)
         cout << ", ";
      else
         cout << "." << endl;
   }

   free_matrix(vindex,1,grid.dim-1);
   for (s = 0; s < 2; s++)
      free_matrix(zindex[s],grid.dim-1,grid.dim-1);
   free_matrix(windex,1,grid.dim-1);
   free_matrix(Du,1,grid.dim-1);
}
// check error
void checkcim345jumpDu(double ***u, double ***S, PBData &pb, GridData &grid)
{
   cout << "[checkcim345jumpDu]" <<endl;
   int i, r, s, t, m, thestatus;
   int tindex[grid.dim], rindex[grid.dim], **zindex = imatrix(grid.dim-1,grid.dim-1);
   double uintp, uintm, Dup[grid.dim], Dum[grid.dim];
   double alpha, thesign, tangent[grid.dim], normal[grid.dim];
   double tmperr, theerr[grid.dim];
   for (t = 0; t < grid.dim; t++)
      theerr[t] = 0.0;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         s = 1;
         rindex[r] = tindex[r]+s;
         if (rindex[r] >= 0 && rindex[r] <= grid.nx[r] && 
             (evalarray(S,tindex) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)
         {
            thestatus = getstatus5(S,tindex,grid);
            if (evalarray(S,tindex) < 0.0)
               thesign = -1.0;
            else
               thesign = 1.0;
            if (thestatus == 2)
               if (thesign < 0.0)
                  if (!globallcim345)
                     getcim3Du(uintm,Dum,tindex,r,s,u,S,pb,grid);
                  else
                     getcim345Du(uintm,Dum,tindex,r,s,u,S,pb,grid);
               else
                  if (!globallcim345)
                     getcim3Du(uintp,Dup,tindex,r,s,u,S,pb,grid);
                  else
                     getcim345Du(uintp,Dup,tindex,r,s,u,S,pb,grid);
            else if (thestatus == 3 || thestatus == 4)
               if (thesign < 0.0)
                  getcim345Du(uintm,Dum,tindex,r,s,u,S,pb,grid);
               else
                  getcim345Du(uintp,Dup,tindex,r,s,u,S,pb,grid);
            else
               cout << "bad for status " << thestatus << " at " 
                    << tindex[0] << " " << tindex[1] << " " << tindex[2] << endl;

            thestatus = getstatus5(S,rindex,grid);
            if (evalarray(S,rindex) < 0.0)
               thesign = -1.0;
            else
               thesign = 1.0;
            if (thestatus == 2)
               if (thesign < 0.0)
                  if (!globallcim345)
                     getcim3Du(uintm,Dum,rindex,r,-s,u,S,pb,grid);
                  else
                     getcim345Du(uintm,Dum,rindex,r,-s,u,S,pb,grid);
               else
                  if (!globallcim345)
                     getcim3Du(uintp,Dup,rindex,r,-s,u,S,pb,grid);
                  else
                     getcim345Du(uintp,Dup,rindex,r,-s,u,S,pb,grid);
            else if (thestatus == 3 || thestatus == 4)
               if (thesign < 0.0)
                  getcim345Du(uintm,Dum,rindex,r,-s,u,S,pb,grid);
               else
                  getcim345Du(uintp,Dup,rindex,r,-s,u,S,pb,grid);
            else
               cout << "bad for status " << thestatus << " at " 
                    << rindex[0] << " " << rindex[1] << " " << rindex[2] << endl;

            getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
            for (t = 0; t < grid.dim; t++)
            {
               tmperr = fabs((Dup[t]-Dum[t])-(getDu(tindex,t,r,s,alpha,1.0,grid)-
                                              getDu(tindex,t,r,s,alpha,-1.0,grid)));
               if (tmperr > theerr[t])
               {
                  theerr[t] = tmperr;
                  for (m = 0; m < grid.dim; m++)
                     zindex[t][m] = tindex[m];
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
   
   cout << "Error in cim345 jump of derivative is " << theerr[0] << " " << theerr[1] 
        << " " << theerr[2] << endl;
   cout << "   at";
   for (t = 0; t < grid.dim; t++)
   {
      for (m = 0; m < grid.dim; m++)
         cout << " " << zindex[t][m];
      cout << " (" << getstatus5(S,zindex[t],grid) << ")";
      if (t < grid.dim-1)
         cout << ", ";
      else
         cout << "." << endl;
   }

   free_matrix(zindex,grid.dim-1,grid.dim-1);
}

void getcim4Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u, 
               double ***S, PBData &pb, GridData &grid)
{
// stencil = 4
   int r, s, t, i, m, n, sk;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];
   double d0[2*grid.dim];
   double ***dcoef[2*grid.dim];
//   double ***D2[grid.dim][grid.dim], D2jumpuxxcoef[grid.dim][grid.dim];
   double ***D2[grid.dim][3], D2jumpuxxcoef[grid.dim][3];
   double **LU, **G;
//   int gamma[grid.dim][2], sk2[grid.dim][grid.dim][4];
   int gamma[grid.dim][2], sk2[grid.dim][3][4];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double alpha, beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[2*grid.dim], **temp2;
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;
   double jumpD1u = 0.0, ***jumpD1ucoef, jumpD1uxcoef[grid.dim], **jumpD1uxxcoef;
   double **jumpD2u, *****jumpD2ucoef, ***jumpD2uxcoef, ***jumpD2uxxcoef;
   double D1u[grid.dim], ***D1ucoef[grid.dim], **D1uxcoef, ***D1uxxcoef;
   double finalalpha, finalD1u[grid.dim], ***finalD1ucoef[grid.dim], **finalD1uxcoef, 
          **finalD1uxxcoef;

   D1uxcoef = matrix(grid.dim-1,grid.dim-1);
   D1uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   finalD1uxcoef = matrix(grid.dim-1,grid.dim-1);
   finalD1uxxcoef = matrix(grid.dim-1,grid.dim-1);
   jumpD1uxxcoef = matrix(grid.dim-1,grid.dim-1);
   temp2 = matrix(grid.dim-1,grid.dim-1);

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

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   G = matrix(2*grid.dim-1,2*grid.dim-1);
   jumpD1ucoef = matrix(4,4,4);
   jumpD2u = matrix(grid.dim-1,grid.dim-1);
   jumpD2uxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   jumpD2uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   jumpD2ucoef = new double ****[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      D1ucoef[r] = matrix(4,4,4);
      finalD1ucoef[r] = matrix(4,4,4);
      jumpD2ucoef[r] = new double ***[grid.dim];
      for (s = 0; s < grid.dim; s++)
      {
         jumpD2ucoef[r][s] = matrix(4,4,4);
         D2[r][s] = matrix(4,4,4);
         D2jumpuxxcoef[r][s] = 0.0;
         temp2[r][s] = 0.0;
      }
   }
   for (r = 0; r < 2*grid.dim; r++)
      dcoef[r] = matrix(4,4,4);

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
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 2;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
            for (s = 0; s < grid.dim; s++)
               D1u[s] = 0.0;
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] < 5)
            {
               for (s = 0; s < grid.dim; s++)
                  setvalarray(D1ucoef[s],sindex,0.0);
               setvalarray(jumpD1ucoef,sindex,0.0);

               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 2;
            
            rindex[r] = index[r]+sk;
// getting interface info
            getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
            beta = 1.0-alpha;
// get sk2 and second derivatives
            for (m = 0; m < grid.dim; m++)
               for (n = m+1; n < grid.dim; n++)
               {
                  if (getsk2(sk2[m][n],m,n,index,S,grid))
                     getD2(D2[m][n],m,n,sk2,sindex,grid);
                  else
                  {
//                     getsk2(sk2[m][n],m,n,rindex,S,grid);
                     if (!getsk2(sk2[m][n],m,n,rindex,S,grid))
                     {
                        cout << "found a non cim4 point" << endl;
                        exit(1);
                     }
                     getD2(D2[m][n],D2jumpuxxcoef[m][n],m,n,r,sk,thesign,sk2,grid);
                  }
               }
            getDu(D1uxcoef,D1uxxcoef,index,r,sk,alpha,thesign,grid);
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha,grid);
            getjumpux(jumpD1u,jumpD1uxcoef,jumpD1uxxcoef,index,r,sk,alpha,thesign,normal,
                      tangent,D1uxcoef,D1uxxcoef,S,pb,grid);
            getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,alpha,
                       thesign,normal,D1uxcoef,D1uxxcoef,D2,D2jumpuxxcoef,S,pb,
                       grid);
            recast(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,jumpD2u,jumpD2ucoef,
                   jumpD2uxcoef,jumpD2uxxcoef,D2,D2jumpuxxcoef,grid);
            if (r == rstar && sk == sstar)
               for (s = 0; s < grid.dim; s++)
                  recast(D1u[s],D1ucoef[s],D1uxcoef[s],D1uxxcoef[s],jumpD2u,jumpD2ucoef,
                         jumpD2uxcoef,jumpD2uxxcoef,D2,D2jumpuxxcoef,grid);
            else
               recast(D1u[r],D1ucoef[r],D1uxcoef[r],D1uxxcoef[r],jumpD2u,jumpD2ucoef,
                      jumpD2uxcoef,jumpD2uxxcoef,D2,D2jumpuxxcoef,grid);
            if (r == rstar && sk == sstar)
            {
               finalalpha = alpha;

               for (s = 0; s < grid.dim; s++)
                  finalD1u[s] = D1u[s];
               for (s = 0; s < grid.dim; s++)
                  for (m = 0; m < grid.dim; m++)
                  {
                     finalD1uxcoef[s][m] = D1uxcoef[s][m];
                     finalD1uxxcoef[s][m] = D1uxxcoef[s][m][m];
                  }
               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 0;
               while (sindex[0] < 5)
               {
                  for (s = 0; s < grid.dim; s++)
                     setvalarray(finalD1ucoef[s],sindex,evalarray(D1ucoef[s],sindex));
                  (sindex[grid.dim-1])++;
                  for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
                  {
                     sindex[i] = 0;
                     (sindex[i-1])++;
                  }
               }
               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 2;
            }
// form d0 and dcoef's rth entry 
            if (globoldwrongcim4)
               d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                                  sk*beta*jumpD1u/grid.dx[r]+
                                                  0.5*beta*beta*jumpD2u[r][r]);
            else
               d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                                  sk*beta*jumpD1u/grid.dx[r]+
                                                  0.5*beta*beta*jumpD2u[r][r])-
                                         sk*D1u[r]/grid.dx[r];

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] < 5)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                           thesign*(sk*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                    0.5*beta*beta*evalarray(jumpD2ucoef[r][r],sindex))-
                           sk*evalarray(D1ucoef[r],sindex)/grid.dx[r]);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 2;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)-
                        1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 2+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)+
                        1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 2;
// form matrix's rth row
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m] = -thesign*
                                           (sk*beta*jumpD1uxxcoef[m][m]/grid.dx[r]+
                                            0.5*beta*beta*jumpD2uxxcoef[r][r][m])+
                                           sk*D1uxxcoef[r][m][m]/grid.dx[r];
            G[grid.dim*(sk+1)/2+r][r] += 0.5*(beta*beta-alpha*alpha);
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m+grid.dim] = -thesign*
                                                    (sk*beta*jumpD1uxcoef[m]/grid.dx[r]+
                                                     0.5*beta*beta*jumpD2uxcoef[r][r][m])+
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
            while (sindex[0] < 5)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,0.0);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 2;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,-1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 2+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = 2;
         }

// solve for uxx and put in matrix
   int j;
   double uxx[grid.dim], ux[grid.dim];
   for (j = 0; j < grid.dim; j++)
   {
      uxx[j] = 0.0;
      ux[j] = 0.0;
   }

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
   while (sindex[0] < 5)
   {
      for (n = 0; n < 2*grid.dim; n++)
         temp[n] = evalarray(dcoef[n],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-2+sindex[s];
      for (j = 0; j < grid.dim; j++)
         uxx[j] += temp[j]*evalarray(u,tindex);
      for (j = 0; j < grid.dim; j++)
//         ux[j] += temp[j+grid.dim]*evalarray(u,tindex);
         if (evalarray(S,tindex) < 0.0)
            ux[j] += temp[j+grid.dim]*getu(tindex,0,0,0.0,-1.0,grid);
         else
            ux[j] += temp[j+grid.dim]*getu(tindex,0,0,0.0,1.0,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 2;

   if (globintorder == 3)
      uint = evalarray(u,index)+sstar*finalalpha*grid.dx[rstar]*ux[rstar]+
             0.5*(finalalpha*grid.dx[rstar])*(finalalpha*grid.dx[rstar])*uxx[rstar];
   else
      uint = evalarray(u,index)+sstar*finalalpha*grid.dx[rstar]*ux[rstar];

   for (s = 0; s < grid.dim; s++)
   {
//      cout << "s = " << s << endl; 
      Du[s] = finalD1u[s];
      for (m = 0; m < grid.dim; m++)
         Du[s] += finalD1uxcoef[s][m]*ux[m]+finalD1uxxcoef[s][m]*uxx[m];
      for (m = 0; m < grid.dim; m++)
         sindex[m] = 0;
      while (sindex[0] < 5)
      {
//         cout << "sindex = " << sindex[0] << " " << sindex[1] << " " << sindex[2] 
//              << endl;
         for (m = 0; m < grid.dim; m++)
            tindex[m] = index[m]-2+sindex[m];
//         cout << "tindex = " << tindex[0] << " " << tindex[1] << " " << tindex[2] 
//              << endl;
         Du[s] += evalarray(finalD1ucoef[s],sindex)*evalarray(u,tindex);
//         cout << "out" << endl;
         (sindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && sindex[i] >= 5; i--)
         {
            sindex[i] = 0;
            (sindex[i-1])++;
         }
      }
      for (m = 0; m < grid.dim; m++)
         sindex[m] = 2;
   }

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

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(G,2*grid.dim-1,2*grid.dim-1);
   free_matrix(jumpD1ucoef,4,4,4);
   free_matrix(jumpD1uxxcoef,grid.dim-1,grid.dim-1);
   free_matrix(jumpD2u,grid.dim-1,grid.dim-1);
   free_matrix(jumpD2uxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(jumpD2uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D1uxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(temp2,grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      free_matrix(D1ucoef[r],4,4,4);
      for (s = 0; s < grid.dim; s++)
      {
         free_matrix(D2[r][s],4,4,4);
         free_matrix(jumpD2ucoef[r][s],4,4,4);
      }
      delete [] jumpD2ucoef[r];
   }
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(dcoef[r],4,4,4);
   delete [] jumpD2ucoef;
}

void getcim5Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u, 
               double ***S, PBData &pb, GridData &grid)
{
   int r, s, t, i, j, m, n, sk, mid = 1, N = 2*mid;
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
   double jumpD1u = 0.0, ***jumpD1ucoef, jumpD1uxcoef[grid.dim], 
          jumpD1uxxcoef[grid.dim];
   double jumpD2u = 0.0, ***jumpD2ucoef, jumpD2uxcoef[grid.dim], 
          jumpD2uxxcoef[grid.dim];
   double ***D1ucoef[grid.dim], **D1uxcoef, **D1uxxcoef, ***D1uxxcoeflarge;
   double finalalpha, ***finalD1ucoef[grid.dim], **finalD1uxcoef, **finalD1uxxcoef;
//   double ***D2[grid.dim][grid.dim], *D2uxcoef[grid.dim][grid.dim], 
//          *D2uxxcoef[grid.dim][grid.dim];
   double ***D2[grid.dim][3], *D2uxcoef[grid.dim][3], 
          *D2uxxcoef[grid.dim][3];
   double ux[grid.dim], uxx[grid.dim];
//   int sk2[grid.dim][grid.dim][4];
   int sk2[grid.dim][3][4];
   

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

   D1uxcoef = matrix(grid.dim-1,grid.dim-1);
   D1uxxcoef = matrix(grid.dim-1,grid.dim-1);
   D1uxxcoeflarge = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   finalD1uxcoef = matrix(grid.dim-1,grid.dim-1);
   finalD1uxxcoef = matrix(grid.dim-1,grid.dim-1);
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   G = matrix(2*grid.dim-1,2*grid.dim-1);
   jumpD1ucoef = matrix(N,N,N);
   jumpD2ucoef = matrix(N,N,N);
   for (r = 0; r < grid.dim; r++)
   {
      D1ucoef[r] = matrix(N,N,N);
      finalD1ucoef[r] = matrix(N,N,N);
      for (s = 0; s < grid.dim; s++)
      {
         D2[r][s] = matrix(N,N,N);
         D2uxcoef[r][s] = new double[grid.dim];
         D2uxxcoef[r][s] = new double[grid.dim];
      }
   }
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
         if (getsk2(sk2[m][n],m,n,index,S,grid) && !globoldcim5)
         {
// make sure sindex is correctly equal to mid
            getD2(D2[m][n],m,n,sk2,sindex,grid);
            for (t = 0; t < grid.dim; t++)
            {
               D2uxcoef[m][n][t] = 0.0;
               D2uxxcoef[m][n][t] = 0.0;
            }
         }
         else
         {
// version with data just used this part of if statement
            getcim5D2(D2[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],m,n,index,mid,S,grid);
         }

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               for (s = 0; s < grid.dim; s++)
                  setvalarray(D1ucoef[s],sindex,0.0);
               setvalarray(jumpD1ucoef,sindex,0.0);

               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
            
            rindex[r] = index[r]+sk;
// getting interface info
            getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
            beta = 1.0-alpha;
// get derivatives
            getDu(D1uxcoef,D1uxxcoeflarge,index,r,sk,alpha,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               recast(D1ucoef[m],D1uxcoef[m],D1uxxcoeflarge[m],D2,D2uxcoef,D2uxxcoef,
                      mid,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = 0; n < grid.dim; n++)
                  if (globoldwrongcim5)
                     D1uxxcoef[m][n] = D1uxxcoeflarge[m][m][n];
                  else
                     D1uxxcoef[m][n] = D1uxxcoeflarge[m][n][n];
            if (r == rstar && sk == sstar)
            {
               finalalpha = alpha;

               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 0;
               while (sindex[0] <= N)
               {
                  for (m = 0; m < grid.dim; m++)
                     setvalarray(finalD1ucoef[m],sindex,evalarray(D1ucoef[m],sindex));

                  (sindex[grid.dim-1])++;
                  for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
                  {
                     sindex[i] = 0;
                     (sindex[i-1])++;
                  }
               }
               for (s = 0; s < grid.dim; s++)
                  sindex[s] = mid;
               for (m = 0; m < grid.dim; m++)
                  for (n = 0; n < grid.dim; n++)
                  {
                     finalD1uxcoef[m][n] = D1uxcoef[m][n];
                     finalD1uxxcoef[m][n] = D1uxxcoef[m][n];
                  }
            }
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha,grid);
            getjumpux(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,index,r,sk,alpha,
                      thesign,normal,tangent,mid,D1ucoef,D1uxcoef,D1uxxcoef,S,pb,grid);
            getjumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,alpha,
                       thesign,normal,mid,D1ucoef,D1uxcoef,D1uxxcoef,D2,D2uxcoef,
                       D2uxxcoef,S,pb,grid);
// form d0 and dcoef's rth entry 
            d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                               sk*beta*jumpD1u/grid.dx[r]+
                                               0.5*beta*beta*jumpD2u);

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

// solve for ux and uxx
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
   {
      Du[s] = 0.0;
      for (m = 0; m < grid.dim; m++)
         sindex[m] = 0;
      while (sindex[0] <= N)
      {
         for (m = 0; m < grid.dim; m++)
            tindex[m] = index[m]-mid+sindex[m];
         Du[s] += evalarray(finalD1ucoef[s],sindex)*evalarray(u,tindex);

         (sindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
         {
            sindex[i] = 0;
            (sindex[i-1])++;
         }
      }
      for (m = 0; m < grid.dim; m++)
         sindex[m] = mid;
      for (m = 0; m < grid.dim; m++)
         Du[s] += finalD1uxcoef[s][m]*ux[m]+finalD1uxxcoef[s][m]*uxx[m];
   }

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(G,2*grid.dim-1,2*grid.dim-1);
   free_matrix(jumpD1ucoef,N,N,N);
   free_matrix(jumpD2ucoef,N,N,N);
   free_matrix(D1uxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1uxxcoeflarge,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(finalD1uxcoef,grid.dim-1,grid.dim-1);
   free_matrix(finalD1uxxcoef,grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      free_matrix(D1ucoef[r],N,N,N);
      free_matrix(finalD1ucoef[r],N,N,N);
      for (s = 0; s < grid.dim; s++)
      {
         free_matrix(D2[r][s],N,N,N);
         delete [] D2uxcoef[r][s];
         delete [] D2uxxcoef[r][s];
      }
   }
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(dcoef[r],N,N,N);
}
// calculate Du- at interface between index and index[rstar]+start
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

void getinterfaceDu(double &uint, double *Du, int *index, int rstar, int sstar,
                    double ***u, double ***S, PBData &pb, GridData &grid)
{
   int thestatus = getstatus5(S,index,grid);

   if (thestatus == 2)
   {
      if (globcim == 2)
//         uint = getinterfacegrad3(Du,u,S,index,rstar,sstar,pb,grid);
         uint = getinterfacegrad4(Du,u,S,index,rstar,sstar,pb,grid);
      else if (!globallcim345)
         getcim3Du(uint,Du,index,rstar,sstar,u,S,pb,grid);
      else
         getcim345Du(uint,Du,index,rstar,sstar,u,S,pb,grid);
   }
   else if (thestatus == 3)
   {
      if (globcim == 2)
//         uint = getinterfacegrad3(Du,u,S,index,rstar,sstar,pb,grid);
         uint = getinterfacegrad4(Du,u,S,index,rstar,sstar,pb,grid);
      else if (globcim == 4)
         getcim4Du(uint,Du,index,rstar,sstar,u,S,pb,grid);
      else if (globcim == 5 || globcim == 345)
         getcim345Du(uint,Du,index,rstar,sstar,u,S,pb,grid);
   }
   else if (thestatus == 4)
   {
      if (globcim == 2)
//         uint = getinterfacegrad3(Du,u,S,index,rstar,sstar,pb,grid);
         uint = getinterfacegrad4(Du,u,S,index,rstar,sstar,pb,grid);
      else if (globcim == 4 || globcim == 5)
         getcim4Du(uint,Du,index,rstar,sstar,u,S,pb,grid);
      else if (globcim == 345)
         getcim345Du(uint,Du,index,rstar,sstar,u,S,pb,grid);
   }
   else if (thestatus == 5)
//      uint = getinterfacegrad3(Du,u,S,index,rstar,sstar,pb,grid);
      uint = getinterfacegrad4(Du,u,S,index,rstar,sstar,pb,grid);
//   getcim345Du(uint,Du,index,rstar,sstar,u,S,pb,grid);
}

void cimexact(SparseElt2**** &A, double ***b, int *index, double ***S, PBData &pb, 
              GridData &grid)
{
   sparse2(index,index,A,1.0,grid);
   if (evalarray(S,index) < 0.0)
      setvalarray(b,index,getu(index,0,0,0.0,-1,grid));
   else
      setvalarray(b,index,getu(index,0,0,0.0,1,grid));
   if (index[0] == 50 && index[1] == 50 && index[2] == 64)
      cout << "exact at this point " << b[index[0]][index[1]][index[2]] << endl;
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

double getexactradius(double thetime, double radius0, double x0, double tol, int Nstep,
                      GridData &grid)
{
   if (globtestnum == 0)
   {
      double x, prevx = x+2.0*(1.0+tol), value, dvalue;
      int step;

      x = x0;
      for (step = 1; step <= Nstep && fabs(x-prevx) > tol; step++)
      {
         value = (log(x)+x*x*(1.0+x*x/4.0))/4.0-
                 (log(radius0)+radius0*radius0*(1.0+radius0*radius0/4.0))/4.0-thetime;//Newtons methodto find root
         dvalue = (1.0/x+x*(2.0+x*x))/4.0;
         prevx = x;
         x -= value/dvalue;
      }
   
      return x;
   }
   else if (globtestnum == 1)
   {
      double epsp = EPSILONP, epsm = EPSILONM;
      return grid.radius0*exp(2.0*(1.0-epsp/epsm)*thetime);
   }
   else if (globtestnum == 2)
   {
      if (thetime == 0.0)
         return radius0;
      else
      {//Heun's Method
         double epsp = EPSILONP, epsm = EPSILONM, temp, vn;
         vn = 2.0*x0*exp(x0*x0)*(1.0-epsp/epsm); //f(r_i)
         temp = x0+grid.dt*vn; // r_i + h f(r_i)
         vn = 2.0*temp*exp(temp*temp)*(1.0-epsp/epsm); // f(r_i + h f(r_i))
         temp = 0.5*(x0+temp)+0.5*grid.dt*vn; //r_i + 0.5 h f(r_i)) + 0.5 h f(r_i + h f(r_i))
         return temp;
      }
   }
   else if (globtestnum == 3 || globtestnum == 4)
   {
      if (thetime == 0.0)
         return radius0;
      else
      {
        //Third-order Strong Stability Preserving Runge-Kutta (SSPRK3), wiki/List_of_Runge-Kutta_methods
         double epsp = EPSILONP, epsm = EPSILONM, tmp, der, rhs, vn; 
         double A = 1.0, B = fabs(1.0-epsp/epsm), C = epsp/epsm;
         vn = x0*(A-C)/sqrt(A*x0*x0+B); //f(r_i) 
         tmp = x0+grid.dt*vn; //x_i + h f(r_i)  
         der = tmp*(A-C)/sqrt(A*tmp*tmp+B); // f(x_i + h f(r_i))
         rhs = vn+der; // f(r_i) + f(x_i + dt f(r_i))
         tmp = x0+0.25*grid.dt*rhs;// r_i +  0.25 h (f(r_i) + f(x_i + dt f(r_i)))
         der = tmp*(A-C)/sqrt(A*tmp*tmp+B); // f(r_i +  0.25 h (f(r_i) + f(x_i + dt f(r_i))))
         rhs = rhs+4.0*der; //  f(r_i) + f(x_i + dt f(r_i)) + 4 f(r_i +  0.25 h (f(r_i) + f(x_i + dt f(r_i))))
         tmp = x0+grid.dt*rhs/6.0; // r_i + 1/6 h (f(r_i) + f(x_i + dt f(r_i)) + 4 f(r_i +  0.25 h (f(r_i) + f(x_i + dt f(r_i))))) 
         return tmp;
      }
   }

   cout << "NOT SUPPOSED TO BE HERE" << endl;
   return 0.0;
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



double getexactresidual(double ***r, SparseElt2 ****A, double ***b, GridData &grid,
                        double ***S, PBData &pb)
{
   int i, m, n, tindex[grid.dim], rindex[grid.dim];
   double bnorm, residual, temp, ehere, thesign;
   SparseElt2 *current;
   double largest = 0.0;
   int lindex[grid.dim];
   
   bnorm = 0.0;
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
            if (evalarray(S,rindex) < 0.0)
               thesign = -1.0;
            else
               thesign = 1.0;
            temp -= 2.0*ehere/(grid.dx[m]*grid.dx[m])*getu(rindex,0,0,0.0,thesign,grid);
            for (n = -1; n <= 1; n += 2)
            {
               rindex[m] = tindex[m]+n;
               if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
               {
                  if (evalarray(S,rindex) < 0.0)
                     thesign = -1.0;
                  else
                     thesign = 1.0;
                  temp -= -ehere/(grid.dx[m]*grid.dx[m])*
                           getu(rindex,0,0,0.0,thesign,grid);
               }
            }
            rindex[m] = tindex[m];
         }
      else
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
         {
            if (evalarray(S,(*current).cindex) < 0.0)
               thesign = -1.0;
            else
               thesign = 1.0;
            temp -= (*current).val*getu((*current).cindex,0,0,0.0,thesign,grid);
         }
      setvalarray(r,tindex,temp);
      if (fabs(temp) > largest)
      {
         largest = fabs(temp);
         for (m = 0; m < grid.dim; m++)
            lindex[m] = tindex[m];
      }
      residual += temp*temp;
      bnorm += evalarray(b,tindex)*evalarray(b,tindex);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);
   residual = sqrt(residual)/bnorm;

   cout << "largest is " << largest << " at " << lindex[0] << " " << lindex[1] << " " 
        << lindex[2] << " with " << getstatus5(S,lindex,grid) << endl;
   return residual;
}

double getexactnormal(int *index, int r, int rstar, int sstar, double alpha, 
                      GridData &grid)
{
   double temp, x[grid.dim];
   int m;

   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];

   return getexactnormalTest(x,r);
}

double getexactDnormal(int *index, int r, int s, int rstar, int sstar, double alpha, 
                       GridData &grid)
{
  cerr<<" exact Dnormal not complete"<<endl;
  exit(1);

   double value, temp, x[grid.dim];
   int m;

   value = -getexactnormal(index,r,rstar,sstar,alpha,grid)*
            getexactnormal(index,s,rstar,sstar,alpha,grid);
   if (r == s)
      value += 1.0;

   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];
   temp = 0.0;
   for (m = 0; m < grid.dim; m++)
      temp += x[m]*x[m];
   temp = sqrt(temp);

   return value/temp;
}
// check error in u with exact solution
void checkexact(double ***S, double radius, PBData &pb, GridData &grid)
{
// use after initializing surface
   int i, r, s, t, tindex[grid.dim], rindex[grid.dim];
   double value, alpha, x[grid.dim], tangent[grid.dim], normal[grid.dim];
   double y[grid.dim];
//   double A = 1.0, B = fabs(1.0-pb.epsilonp/pb.epsilonm),
//          C = pb.epsilonp/pb.epsilonm, D = radius*radius*(A-C)+B;
   double theerr[4], thesign, ehere, tau, sigma;

   for (i = 0; i < grid.dim; i++)
      theerr[i] = 0.0;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (evalarray(S,tindex) < 0.0)
      {
         thesign = -1.0;
         ehere = pb.epsilonm;
      }
      else
      {
         thesign = 1.0;
         ehere = pb.epsilonp;
      }
      value = 0.0;
      for (r = 0; r < grid.dim; r++)
         value += getD2u(tindex,r,r,0,0,0.0,thesign,grid);
      value *= ehere;
      value += getf(tindex,0,0,0.0,thesign,pb,grid);
      if (fabs(value) > theerr[0])
         theerr[0] = fabs(value);
     
      // eta2d << getu(tindex,0,0,0.0,thesign,grid) << " ";

      s = 1;
      for (r = 0; r < grid.dim; r++)
      {
         for (i = 0; i < grid.dim; i++)
            rindex[i] = tindex[i];
         rindex[r] = tindex[r]+s;//look at tindex + (0,0,1) 
         if (rindex[r] <= grid.nx[r] && 
             (evalarray(S,tindex) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)//exactly one of tindex or rindex in inside
         {
            sub2coord(x,tindex,grid);
            getinterfaceinfo(alpha,tangent,normal,S,tindex,rindex,grid);
            x[r] += s*alpha*grid.dx[r];

//            gettau(tau,tindex,r,s,alpha,grid);
            tau = 0.0;
            value = getu(tindex,r,s,alpha,1.0,grid)-
                    getu(tindex,r,s,alpha,-1.0,grid)-tau;
            if (fabs(value) > theerr[1])
               theerr[1] = fabs(value);

            sigma = 0.0;
            value = 0.0;
            for (t = 0; t < grid.dim; t++)
               value += (pb.epsilonp*getDu(tindex,t,r,s,alpha,1.0,grid)-
                         pb.epsilonm*getDu(tindex,t,r,s,alpha,-1.0,grid))*
                        getexactnormal(tindex,t,r,s,alpha,grid);
            value -= sigma;
            if (fabs(value) > theerr[2])
               theerr[2] = fabs(value);
         }
      }

      for (r = 0; r < grid.dim && tindex[r] > 0 && tindex[r] < grid.nx[r]; r++); //if r==3, then tindex out of bound
      if (r < grid.dim)//if tindex inbound
      {
         sub2coord(x,tindex,grid);
         value = getu(tindex,0,0,0.0,thesign,grid)-DBC(x,grid.dim,0.0);
         if (fabs(value) > theerr[3])
            theerr[3] = fabs(value);
      }

      // (tindex[grid.dim-1])++;
      // for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      // {
      //    tindex[i] = 0;
      //    (tindex[i-1])++;
      //    if (i == grid.dim-1)
      //       eta2d << endl;
      // }
   }

   cout << "Max error in interior: " << theerr[0] << endl;
   cout << "Max errors at boundary: " << theerr[1] << " and " << theerr[2] << endl;
   cout << "Max errors at computational boundary: " << theerr[3] << endl;
}

void getinterfaceinfo(double *normal, double *tangent1, double *tangent2, double **Dn,
                      double &tau, double *Dtau, double **D2tau, double &sigma, 
                      double *Dsigma, double &jumpfe, double *x, double ***S,
                      PBData &pb, GridData &grid)
{
   int i, r, s, t, m, n; 
   int index[grid.dim], tindex[grid.dim], rindex[grid.dim], sindex[grid.dim];
   double grad, length, tol = 1.0e-14, val, tempx[grid.dim];
   double Dn2[grid.dim][grid.dim];
   double ***temp = matrix(1,1,1);
   GridData tempgrid;

   for (r = 0; r < grid.dim; r++)
   {
      index[r] = (x[r]-grid.a[r])/grid.dx[r];
      if (index[r] >= grid.nx[r])
         index[r] = grid.nx[r]-1;
      else if (index[r] < 0)
         index[r] = 0;
      tempgrid.nx[r] = 1;
      tempgrid.a[r] = 0.0;
      tempgrid.dx[r] = grid.dx[r];
   }
   sub2coord(tempx,index,grid);
   for (r = 0; r < grid.dim; r++)
      tempx[r] = x[r]-tempx[r];
   tempgrid.tol = grid.tol;
   for (r = 0; r < grid.dim; r++)
   {
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= 1)
      {
         for (t = 0; t < grid.dim; t++)
         {
            sindex[t] = index[t]+tindex[t];
            rindex[t] = sindex[t];
         }
         val = 0.0;
         for (s = -1; s <= 1; s += 2) 
         {
            rindex[r] = min(max(sindex[r]+s,0),grid.nx[r]);
            val += s*evalarray(S,rindex);
         }
         val /= 2.0*grid.dx[r];
         setvalarray(temp,tindex,val);
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > 1; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      normal[r] = bilinearinterp(tempx,temp,tempgrid);
   }
   grad = sqrt(getdotprod(normal,normal,grid.dim));
   for (r = 0; r < grid.dim; r++)
      normal[r] /= grad;

// use exact here
//   for (r = 0; r < grid.dim; r++)
//      normal[r] = x[r]/sqrt(getdotprod(x,x,grid.dim));

   length = 0.0;
   t = -1;
   for (r = 0; r < grid.dim; r++)
   {   
      for (s = 0; s < grid.dim; s++)
         tangent2[s] = 0.0;
      tangent2[r] = 1.0;
      project(tangent2,normal,tangent2,grid.dim);
      if (sqrt(getdotprod(tangent2,tangent2,grid.dim)) > length)
      {
         length = sqrt(getdotprod(tangent2,tangent2,grid.dim));
         for (s = 0; s < grid.dim; s++)
            tangent1[s] = tangent2[s]/length;
      }
   }
   getunitcrossprod(tangent2,normal,tangent1);

   for (r = 0; r < grid.dim; r++)
      for (m = r; m < grid.dim; m++)
      {
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= 1)
         {
            for (t = 0; t < grid.dim; t++)
            {
               sindex[t] = index[t]+tindex[t];
               rindex[t] = sindex[t];
            }
            val = 0.0;
            if (r == m)
            {
               for (s = -1; s <= 1; s += 2) 
               {
                  rindex[r] = min(max(sindex[r]+s,0),grid.nx[r]);
                  val += evalarray(S,rindex)-evalarray(S,sindex);
               }
               val /= grid.dx[r]*grid.dx[r];
            }
            else
            {
               for (s = -1; s <= 1; s += 2)
               {
                  rindex[r] = min(max(sindex[r]+s,0),grid.nx[r]);
                  for (t = -1; t <= 1; t += 2)
                  {
                     rindex[m] = min(max(sindex[m]+t,0),grid.nx[r]);
                     val += s*t*evalarray(S,rindex);
                  }
               }
               val /= 4.0*grid.dx[r]*grid.dx[m];
            }
            setvalarray(temp,tindex,val);
   
            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > 1; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }
         Dn2[r][m] = bilinearinterp(tempx,temp,tempgrid);
         Dn2[m][r] = Dn2[r][m];
      }

   for (r = 0; r < grid.dim; r++)
      for (m = 0; m < grid.dim; m++)
      {
         Dn[r][m] = Dn2[r][m]/grad;
         for (n = 0; n < grid.dim; n++)
            Dn[r][m] -= normal[r]*Dn2[m][n]*normal[n]/grad;
      }

// use exact Dn
//   for (r = 0; r < grid.dim; r++)
//      for (s = 0; s < grid.dim; s++)
//      {
//         Dn[r][s] = -x[r]*x[s]/sqrt(getdotprod(x,x,grid.dim))/getdotprod(x,x,grid.dim);
//         if (r == s)
//            Dn[r][s] += 1.0/sqrt(getdotprod(x,x,grid.dim));
//      }

   jumpfe = getf(x,1,pb,grid)/pb.epsilonp-getf(x,-1,pb,grid)/pb.epsilonm;
   tau = gettau(x,pb,grid);
   getDtau(Dtau,x,pb,grid);
   getD2tau(D2tau,x,pb,grid);
   sigma = getsigma(x,normal,pb,grid);
   getDsigma(Dsigma,x,normal,Dn,pb,grid);

/*
   double tempval;
   cout << "testing" << endl;
   for (r = 0; r < grid.dim; r++)
      for (s = 0; s < grid.dim; s++)
      {
         tempval = -x[r]*x[s]/sqrt(getdotprod(x,x,grid.dim))/getdotprod(x,x,grid.dim);
         if (r == s)
            tempval += 1.0/sqrt(getdotprod(x,x,grid.dim));
         cout << tempval << " " << Dn[r][s] << endl;
      }
*/

   free_matrix(temp,1,1,1);
}

void getiimjumps(double &up0, double &upum, double *uxp0, double **uxpuxm, 
                 double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm,
                 int *index, int rstar, int sk, double alpha, int thesign, 
                 double *normal, double ***S, PBData &pb, GridData &grid)
{
   int r, s, n, m;
   double **LU, **B, **Dn; 
   double Dndott[grid.dim], Dndots[grid.dim], dotDndot[grid.dim];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double tangent1[grid.dim], tangent2[grid.dim];
   double ethere, ehere, jumpfe;
   double temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim];
   double Dtau[grid.dim], **D2tau, dotD2taudot[grid.dim], Dtaudot[grid.dim];
   int two2one[grid.dim][grid.dim];
//   double aehere, aethere, jumpae;

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   B = matrix(2*grid.dim-1,2*grid.dim-1);
   Dn = matrix(grid.dim-1,grid.dim-1);
   D2tau = matrix(grid.dim-1,grid.dim-1);

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
//   jumpae = thesign*(aehere-aethere);
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

   up0 = -thesign*tau;
   upum = 1.0;

   for (r = 0; r < grid.dim; r++)
      uxp0[r] = -thesign*(normal[r]*sigma/ethere+tangent1[r]*Dtaudot[1]+
                                                 tangent2[r]*Dtaudot[2]);
   for (r = 0; r < grid.dim; r++)
      for (s = 0; s < grid.dim; s++)
         uxpuxm[r][s] = normal[r]*normal[s]*(1.0+(ehere-ethere)/ethere)+
                        tangent1[r]*tangent1[s]+tangent2[r]*tangent2[s];

   gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);

   for (r = 0; r < grid.dim; r++)
      temp[r] = -thesign*(dotD2taudot[r]+(sigma/ethere-Dtaudot[0])*dotDndot[r]);
   temp[grid.dim] = -thesign*(getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                              Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2]);
   temp[grid.dim+1] = -thesign*(getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                                Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1]);
   temp[grid.dim+2] = thesign*jumpfe;
   forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
   for (r = 0; r < grid.dim; r++)
      for (s = r; s < grid.dim; s++)
         uxxp0[r][s] = temp[two2one[r][s]];

   for (m = 0; m < grid.dim; m++)
   {
      for (r = 0; r < grid.dim; r++)
         temp[r] = (ehere-ethere)/ethere*dotDndot[r]*normal[m];
      temp[grid.dim] = (ehere-ethere)/ethere*Dndots[m];
      temp[grid.dim+1] = (ehere-ethere)/ethere*Dndott[m];
      temp[grid.dim+2] = 0.0;
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      for (r = 0; r < grid.dim; r++)
         for (s = r; s < grid.dim; s++)
            uxxpuxm[r][s][m] = temp[two2one[r][s]];
   }

   for (m = 0; m < grid.dim; m++)
      for (n = m; n < grid.dim; n++)
      {
         for (r = 0; r < grid.dim; r++)
            temp[r] = B[r][two2one[m][n]];
         for (r = grid.dim; r < grid.dim+2; r++)
            temp[r] = (1.0+(ehere-ethere)/ethere)*B[r][two2one[m][n]];
         if (m != n)
            temp[grid.dim+2] = 0.0;
         else
            temp[grid.dim+2] = 1.0+(ehere-ethere)/ethere;
         forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
         for (r = 0; r < grid.dim; r++)
            for (s = r; s < grid.dim; s++)
               uxxpuxxm[r][s][m][n] = temp[two2one[r][s]];
      }

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
}

void getiimjumps(double &up0, double &upum, double *uxp0, double **uxpuxm, 
                 double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm,
                 double *x, int thesign, double ***S, PBData &pb, GridData &grid)
{
   int r, s, n, m;
   double **LU, **B, **Dn; 
   double Dndott[grid.dim], Dndots[grid.dim], dotDndot[grid.dim];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double tangent1[grid.dim], tangent2[grid.dim];
   double ethere, ehere, jumpfe;
   double temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim];
   double Dtau[grid.dim], **D2tau, dotD2taudot[grid.dim], Dtaudot[grid.dim], 
          normal[grid.dim];
   int two2one[grid.dim][grid.dim];
//   double aehere, aethere, jumpae;

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   B = matrix(2*grid.dim-1,2*grid.dim-1);
   Dn = matrix(grid.dim-1,grid.dim-1);
   D2tau = matrix(grid.dim-1,grid.dim-1);

   if (thesign < 0)
   {
      ehere = pb.epsilonm;
      ethere = pb.epsilonp;
   }
   else
   {
      ehere = pb.epsilonp;
      ethere = pb.epsilonm;
   }

   getinterfaceinfo(normal,tangent1,tangent2,Dn,tau,Dtau,D2tau,sigma,Dsigma,jumpfe,x,S,
                    pb,grid);
// form Dn dot product with various vectors
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

   up0 = -thesign*tau;
   upum = 1.0;

//   cout << "Check up to um: "
//        << fabs(getu(x,-thesign,grid)-(upum*getu(x,thesign,grid)+up0)) << endl;

   for (r = 0; r < grid.dim; r++)
      uxp0[r] = -thesign*(normal[r]*sigma/ethere+tangent1[r]*Dtaudot[1]+
                                                 tangent2[r]*Dtaudot[2]);
   for (r = 0; r < grid.dim; r++)
      for (s = 0; s < grid.dim; s++)
         uxpuxm[r][s] = normal[r]*normal[s]*(1.0+(ehere-ethere)/ethere)+
                        tangent1[r]*tangent1[s]+tangent2[r]*tangent2[s];

//   cout << "Check uxp to uxm: "
//        << fabs(getDu(x,0,-thesign,grid)-(uxpuxm*getDu(x,0,thesign,grid)+uxp0)) << " "
//        << fabs(getDu(x,1,-thesign,grid)-(uxpuxm*getDu(x,1,thesign,grid)+uxp0)) << " "
//        << fabs(getDu(x,2,-thesign,grid)-(uxpuxm*getDu(x,2,thesign,grid)+uxp0)) << endl;

   gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);

   for (r = 0; r < grid.dim; r++)
      temp[r] = dotD2taudot[r]+(sigma/ethere-Dtaudot[0])*dotDndot[r];
   temp[grid.dim] = getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                    Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2];
   temp[grid.dim+1] = getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                      Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1];
   temp[grid.dim+2] = -jumpfe;
   forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
   for (r = 0; r < grid.dim; r++)
      for (s = r; s < grid.dim; s++)
         uxxp0[r][s] = temp[two2one[r][s]];

   for (m = 0; m < grid.dim; m++)
   {
      for (r = 0; r < grid.dim; r++)
         temp[r] = thesign*(ethere-ehere)/ethere*dotDndot[r]*normal[m];
      temp[grid.dim] = thesign*(ethere-ehere)/ethere*Dndott[m];
      temp[grid.dim+1] = thesign*(ethere-ehere)/ethere*Dndots[m];
      temp[grid.dim+2] = 0.0;
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      for (r = 0; r < grid.dim; r++)
         for (s = r; s < grid.dim; s++)
            uxxpuxm[r][s][m] = temp[two2one[r][s]];
   }

   for (m = 0; m < grid.dim; m++)
      for (n = m; n < grid.dim; n++)
      {
         for (r = 0; r < grid.dim; r++)
            temp[r] = 0.0;
         for (r = grid.dim; r < grid.dim+2; r++)
            temp[r] = thesign*(ethere-ehere)/ethere*B[r][two2one[m][n]];
         temp[grid.dim+2] = 0.0;
         forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
         for (r = 0; r < grid.dim; r++)
            for (s = r; s < grid.dim; s++)
               uxxpuxxm[r][s][m][n] = temp[two2one[r][s]];
      }

   for (r = 0; r < grid.dim; r++)
      for (s = r; s < grid.dim; s++)
      {
         uxxp0[r][s] *= -thesign;
         for (m = 0; m < grid.dim; m++)
            uxxpuxm[r][s][m] *= -thesign;
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               uxxpuxxm[r][s][m][n] *= -thesign;
         uxxpuxxm[r][s][r][s] += 1.0;
      }

   if (globdebug)
   {
      cout << "in iimjumps" << endl;
      double tempval;
      cout << getu(x,-thesign,grid)-(upum*getu(x,thesign,grid)+up0) << endl;
      for (r = 0; r < grid.dim; r++)
      {
         tempval = uxp0[r];
         for (s = 0; s < grid.dim; s++)
            tempval += uxpuxm[r][s]*getDu(x,s,thesign,grid);
         cout << getDu(x,r,-thesign,grid)-tempval << endl;
      }
      for (r = 0; r < grid.dim; r++)
         for (s = r; s < grid.dim; s++)
         {
            tempval = uxxp0[r][s];
            for (m = 0; m < grid.dim; m++)
               tempval += uxxpuxm[r][s][m]*getDu(x,m,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  tempval += uxxpuxxm[r][s][m][n]*getD2u(x,m,n,thesign,grid);
            cout << getD2u(x,r,s,-thesign,grid)-tempval << endl;
//         cout << "   " 
//              << getD2u(x,r,s,1.0,grid)-getD2u(x,r,s,-1.0,grid)-tempval << endl;
//         cout << "   " 
//              << -getD2u(x,r,s,1.0,grid)+getD2u(x,r,s,-1.0,grid)-tempval << endl;
         }
   }

/*
   double tempvec[2*grid.dim];
   tempvec[0] = (sigma/ethere+thesign*(ethere-ehere)/ethere*
                 (getDu(x,0,thesign,grid)*normal[0]+getDu(x,1,thesign,grid)*normal[1]+
                  getDu(x,2,thesign,grid)*normal[2]))*dotDndot[0];
   tempvec[1] = (sigma/ethere+thesign*(ethere-ehere)/ethere*
                 (getDu(x,0,thesign,grid)*normal[0]+getDu(x,1,thesign,grid)*normal[1]+
                  getDu(x,2,thesign,grid)*normal[2]))*dotDndot[1];
   tempvec[2] = (sigma/ethere+thesign*(ethere-ehere)/ethere*
                 (getDu(x,0,thesign,grid)*normal[0]+getDu(x,1,thesign,grid)*normal[1]+
                  getDu(x,2,thesign,grid)*normal[2]))*dotDndot[2];
   tempvec[3] = getdotprod(Dsigma,tangent1,grid.dim)/ethere+
                thesign*(ethere-ehere)/ethere*
                (getDu(x,0,thesign,grid)*Dndott[0]+getDu(x,1,thesign,grid)*Dndott[1]+
                  getDu(x,2,thesign,grid)*Dndott[2])+
                thesign*(ethere-ehere)/ethere*
                (normal[0]*getD2u(x,0,0,thesign,grid)*tangent1[0]+
                 normal[1]*getD2u(x,1,0,thesign,grid)*tangent1[0]+
                 normal[2]*getD2u(x,2,0,thesign,grid)*tangent1[0]+
                 normal[0]*getD2u(x,0,1,thesign,grid)*tangent1[1]+
                 normal[1]*getD2u(x,1,1,thesign,grid)*tangent1[1]+
                 normal[2]*getD2u(x,2,1,thesign,grid)*tangent1[1]+
                 normal[0]*getD2u(x,0,2,thesign,grid)*tangent1[2]+
                 normal[1]*getD2u(x,1,2,thesign,grid)*tangent1[2]+
                 normal[2]*getD2u(x,2,2,thesign,grid)*tangent1[2]);
   tempvec[4] = getdotprod(Dsigma,tangent2,grid.dim)/ethere+
                thesign*(ethere-ehere)/ethere*
                (getDu(x,0,thesign,grid)*Dndots[0]+getDu(x,1,thesign,grid)*Dndots[1]+
                  getDu(x,2,thesign,grid)*Dndots[2])+
                thesign*(ethere-ehere)/ethere*
                (normal[0]*getD2u(x,0,0,thesign,grid)*tangent2[0]+
                 normal[1]*getD2u(x,1,0,thesign,grid)*tangent2[0]+
                 normal[2]*getD2u(x,2,0,thesign,grid)*tangent2[0]+
                 normal[0]*getD2u(x,0,1,thesign,grid)*tangent2[1]+
                 normal[1]*getD2u(x,1,1,thesign,grid)*tangent2[1]+
                 normal[2]*getD2u(x,2,1,thesign,grid)*tangent2[1]+
                 normal[0]*getD2u(x,0,2,thesign,grid)*tangent2[2]+
                 normal[1]*getD2u(x,1,2,thesign,grid)*tangent2[2]+
                 normal[2]*getD2u(x,2,2,thesign,grid)*tangent2[2]);
   tempvec[5] = -(getf(x,1,pb,grid)/pb.epsilonp-getf(x,-1,pb,grid)/pb.epsilonm);
   forwardbacksub0(tempvec,tempvec,LU,PLR,PLC,2*grid.dim-1);
   cout << "testing" << endl;
   cout << tempvec[0] << " " << getD2u(x,0,0,1.0,grid)-
                                getD2u(x,0,0,-1.0,grid) << " " 
        << tempvec[0]-(getD2u(x,0,0,1.0,grid)-getD2u(x,0,0,-1.0,grid)) << endl;
   cout << tempvec[1] << " " << getD2u(x,1,1,1.0,grid)-
                                getD2u(x,1,1,-1.0,grid) << " "
        << tempvec[1]-(getD2u(x,1,1,1.0,grid)-getD2u(x,1,1,-1.0,grid)) << endl;
   cout << tempvec[2] << " " << getD2u(x,2,2,1.0,grid)-
                                getD2u(x,2,2,-1.0,grid) << " "
        << tempvec[2]-(getD2u(x,2,2,1.0,grid)-getD2u(x,2,2,-1.0,grid)) << endl;
   cout << tempvec[3] << " " << getD2u(x,0,1,1.0,grid)-
                                getD2u(x,0,1,-1.0,grid) << " "
        << tempvec[3]-(getD2u(x,0,1,1.0,grid)-getD2u(x,0,1,-1.0,grid)) << endl;
   cout << tempvec[4] << " " << getD2u(x,1,2,1.0,grid)-
                                getD2u(x,1,2,-1.0,grid) << " "
        << tempvec[4]-(getD2u(x,1,2,1.0,grid)-getD2u(x,1,2,-1.0,grid)) << endl;
   cout << tempvec[5] << " " << getD2u(x,0,2,1.0,grid)-
                                getD2u(x,0,2,-1.0,grid) << " "
        << tempvec[5]-(getD2u(x,0,2,1.0,grid)-getD2u(x,0,2,-1.0,grid)) << endl;
*/

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
}

/*
void getiimstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                     double up0, double upum, double *uxp0, double **uxpuxm, 
                     double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm, 
                     double ***S, PBData &pb, GridData &grid)
{

   double **z = matrix(N-1,grid.dim-1);
   int *signs = new int[N];
   int two2one[grid.dim][grid.dim];
   int i, r, s, m, n;

   for (i = 0; i < N; i++)
   {
      sub2coord(z[i],index[i],grid);
      for (r = 0; r < grid.dim; r++)
         z[i][r] = z[i][r]-x[r];
      if (evalarray(S,index[i]) < 0.0)
         signs[i] = -1;
      else
         signs[i] = 1;
   }
   for (n = 0; n < grid.dim; n++)
      two2one[n][n] = n;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (m = 0; m < grid.dim; m++)
      for (n = 0; n < m; n++)
         two2one[m][n] = two2one[n][m];

   for (i = 0; i < N; i++)
   {
      if (signs[i] == thesign)
      {
         A[0][i] = 0.0;
         A[1][i] = 1.0;
         for (m = 0; m < grid.dim; m++)
            A[2+m][i] = z[i][m];
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               if (m == n)
                  A[2+grid.dim+two2one[m][n]][i] = 0.5*z[i][m]*z[i][n];
               else
                  A[2+grid.dim+two2one[m][n]][i] = z[i][m]*z[i][n];
      }
      else
      {
         A[0][i] = up0;
         for (r = 0; r < grid.dim; r++)
            A[0][i] += z[i][r]*uxp0[r];
         for (r = 0; r < grid.dim; r++)
            for (s = r; s < grid.dim; s++)
               if (r == s)
                  A[0][i] += 0.5*z[i][r]*z[i][s]*uxxp0[r][s];
               else
                  A[0][i] += z[i][r]*z[i][s]*uxxp0[r][s];

         A[1][i] = upum;

         for (m = 0; m < grid.dim; m++)
         {
            A[2+m][i] = 0.0;
            for (r = 0; r < grid.dim; r++)
               A[2+m][i] += z[i][r]*uxpuxm[m][r];
            for (r = 0; r < grid.dim; r++)
               for (s = r; s < grid.dim; s++)
                  if (r == s)
                     A[2+m][i] += 0.5*z[i][r]*z[i][s]*uxxpuxm[r][s][m];
                  else
                     A[2+m][i] += z[i][r]*z[i][s]*uxxpuxm[r][s][m];
         }

         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
            {
               A[2+grid.dim+two2one[m][n]][i] = 0.0;
               for (r = 0; r < grid.dim; r++)
                  for (s = r; s < grid.dim; s++)
                     if (r == s)
                        A[2+grid.dim+two2one[m][n]][i] += 0.5*z[i][r]*z[i][s]*
                                                          uxxpuxxm[r][s][m][n];
                     else
                        A[2+grid.dim+two2one[m][n]][i] = z[i][r]*z[i][s]*
                                                         uxxpuxxm[r][s][m][n];
            }
      }
   }
//   gecp0(LU,PLR,PLC,A,M-1,N-1);
//   for (r = 0; r < grid.dim; r++)
//   {
//      b[2+grid.dim+two2one[r][r]] = 1.0;
//      gecp0(LU,PLR,PLC,A,M-1,N-1);
//      forwardbacksub0(c[r],b,LU,PLR,PLC,M-1,N-1);
//      b[2+grid.dim+two2one[r][r]] = 0.0;
//   }

   free_matrix(z,N-1,grid.dim-1);
   delete [] signs;
}
*/

void getiimstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                     double up0, double upum, double *uxp0, double **uxpuxm, 
                     double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm, 
                     double ***S, PBData &pb, GridData &grid)
{

   double **z = matrix(N-1,grid.dim-1);
   int *signs = new int[N];
   int two2one[grid.dim][grid.dim];
   int i, r, s, m, n;

   cout << "main point " << x[0] << " " << x[1] << " " << x[2] << " has sign " 
        << thesign << endl;
   for (i = 0; i < N; i++)
   {
      sub2coord(z[i],index[i],grid);
      for (r = 0; r < grid.dim; r++)
         z[i][r] = z[i][r]-x[r];
      if (evalarray(S,index[i]) < 0.0)
         signs[i] = -1;
      else
         signs[i] = 1;
      cout << "point at " << index[i][0] << " " << index[i][1] << " " << index[i][2] 
           << " has sign " << signs[i] << " and dist " 
           << sqrt(getdotprod(z[i],z[i],grid.dim)) << endl;
   }
   for (n = 0; n < grid.dim; n++)
      two2one[n][n] = n;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (m = 0; m < grid.dim; m++)
      for (n = 0; n < m; n++)
         two2one[m][n] = two2one[n][m];

   double tempval;

   for (i = 0; i < N; i++)
   {
      if (signs[i] == thesign)
      {
         A[0][i] = 0.0;
         A[1][i] = 1.0;
         for (m = 0; m < grid.dim; m++)
            A[2+m][i] = z[i][m];
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               if (m == n)
                  A[2+grid.dim+two2one[m][n]][i] = 0.5*z[i][m]*z[i][n];
               else
                  A[2+grid.dim+two2one[m][n]][i] = z[i][m]*z[i][n];
         tempval = A[0][i];
         tempval += A[1][i]*getu(x,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            tempval += A[2+m][i]*getDu(x,m,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               tempval += A[2+grid.dim+two2one[m][n]][i]*getD2u(x,m,n,thesign,grid);
         cout << "Taylor series diff = " << getu(index[i],0,0,0.0,thesign,grid)-tempval 
              << endl;
      }
      else
      {
         A[0][i] = up0;
         for (r = 0; r < grid.dim; r++)
            A[0][i] += z[i][r]*uxp0[r];
         for (r = 0; r < grid.dim; r++)
            for (s = r; s < grid.dim; s++)
               if (r == s)
                  A[0][i] += 0.5*z[i][r]*z[i][s]*uxxp0[r][s];
               else
                  A[0][i] += z[i][r]*z[i][s]*uxxp0[r][s];

         A[1][i] = upum;

         for (m = 0; m < grid.dim; m++)
         {
            A[2+m][i] = 0.0;
            for (r = 0; r < grid.dim; r++)
               A[2+m][i] += z[i][r]*uxpuxm[r][m];
            for (r = 0; r < grid.dim; r++)
               for (s = r; s < grid.dim; s++)
                  if (r == s)
                     A[2+m][i] += 0.5*z[i][r]*z[i][s]*uxxpuxm[r][s][m];
                  else
                     A[2+m][i] += z[i][r]*z[i][s]*uxxpuxm[r][s][m];
         }

         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
            {
               A[2+grid.dim+two2one[m][n]][i] = 0.0;
               for (r = 0; r < grid.dim; r++)
                  for (s = r; s < grid.dim; s++)
                     if (r == s)
                        A[2+grid.dim+two2one[m][n]][i] += 0.5*z[i][r]*z[i][s]*
                                                          uxxpuxxm[r][s][m][n];
                     else
                        A[2+grid.dim+two2one[m][n]][i] += z[i][r]*z[i][s]*
                                                          uxxpuxxm[r][s][m][n];
            }

         tempval = A[0][i];
         tempval += A[1][i]*getu(x,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            tempval += A[2+m][i]*getDu(x,m,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               tempval += A[2+grid.dim+two2one[m][n]][i]*getD2u(x,m,n,thesign,grid);
         cout << "Taylor series diff = " << getu(index[i],0,0,0.0,-thesign,grid)-tempval 
              << endl;

         tempval = 0.0;
         tempval += getu(x,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            tempval += z[i][m]*getDu(x,m,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               if (m == n)
                  tempval += 0.5*z[i][m]*z[i][n]*getD2u(x,m,n,thesign,grid);
               else
                  tempval += z[i][m]*z[i][n]*getD2u(x,m,n,thesign,grid);
         cout << "   on other side = " << getu(index[i],0,0,0.0,thesign,grid)-tempval 
              << endl;
      }
   }

   cout << "the " << M << "x" << N << " matrix = " << endl;
   for (m = 0; m < M; m++)
   {
      for (n = 0; n < N; n++)
         cout << A[m][n] << " ";
      cout << endl;
   }

   free_matrix(z,N-1,grid.dim-1);
   delete [] signs;
}

void getiimstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                     char yesC, double up0, double upum, double *uxp0, 
                     double **uxpuxm, double **uxxp0, double ***uxxpuxm, 
                     double ****uxxpuxxm, double ***S, PBData &pb, GridData &grid)
{

   double **z = matrix(N-1,grid.dim-1);
   int *signs = new int[N];
   int two2one[grid.dim][grid.dim];
   int i, r, s, m, n;

   if (globdebug)
      cout << "main point " << x[0] << " " << x[1] << " " << x[2] << " has sign " 
           << thesign << endl;
   for (i = 0; i < N; i++)
   {
      sub2coord(z[i],index[i],grid);
      for (r = 0; r < grid.dim; r++)
         z[i][r] = z[i][r]-x[r];
      if (evalarray(S,index[i]) < 0.0)
         signs[i] = -1;
      else
         signs[i] = 1;
      if (globdebug)
         cout << "point at " << index[i][0] << " " << index[i][1] << " " << index[i][2] 
              << " has sign " << signs[i] << " and dist " 
              << sqrt(getdotprod(z[i],z[i],grid.dim)) << endl;
   }
   for (n = 0; n < grid.dim; n++)
      two2one[n][n] = n;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (m = 0; m < grid.dim; m++)
      for (n = 0; n < m; n++)
         two2one[m][n] = two2one[n][m];

   double tempval, temperr[10];

   if (yesC)
   {
      A[0][N] = 1.0;
      for (i = 1; i < M; i++)
         A[i][N] = 0.0;
   }
   for (i = 0; i < N; i++)
   {
      if (signs[i] == thesign)
      {
         A[0][i] = 0.0;
         A[1][i] = 1.0;
         for (m = 0; m < grid.dim; m++)
            A[2+m][i] = z[i][m];
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               if (m == n)
                  A[2+grid.dim+two2one[m][n]][i] = 0.5*z[i][m]*z[i][n];
               else
                  A[2+grid.dim+two2one[m][n]][i] = z[i][m]*z[i][n];
         if (globdebug)
         {
            tempval = A[0][i];
            tempval += A[1][i]*getu(x,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               tempval += A[2+m][i]*getDu(x,m,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  tempval += A[2+grid.dim+two2one[m][n]][i]*getD2u(x,m,n,thesign,grid);
            cout << "Taylor series diff = " 
                 << getu(index[i],0,0,0.0,thesign,grid)-tempval << endl;
         }
      }
      else
      {
         A[0][i] = up0;
         for (r = 0; r < grid.dim; r++)
            A[0][i] += z[i][r]*uxp0[r];
         for (r = 0; r < grid.dim; r++)
            for (s = r; s < grid.dim; s++)
               if (r == s)
                  A[0][i] += 0.5*z[i][r]*z[i][s]*uxxp0[r][s];
               else
                  A[0][i] += z[i][r]*z[i][s]*uxxp0[r][s];

         A[1][i] = upum;

         for (m = 0; m < grid.dim; m++)
         {
            A[2+m][i] = 0.0;
            for (r = 0; r < grid.dim; r++)
               A[2+m][i] += z[i][r]*uxpuxm[r][m];
            for (r = 0; r < grid.dim; r++)
               for (s = r; s < grid.dim; s++)
                  if (r == s)
                     A[2+m][i] += 0.5*z[i][r]*z[i][s]*uxxpuxm[r][s][m];
                  else
                     A[2+m][i] += z[i][r]*z[i][s]*uxxpuxm[r][s][m];
         }

         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
            {
               A[2+grid.dim+two2one[m][n]][i] = 0.0;
               for (r = 0; r < grid.dim; r++)
                  for (s = r; s < grid.dim; s++)
                     if (r == s)
                        A[2+grid.dim+two2one[m][n]][i] += 0.5*z[i][r]*z[i][s]*
                                                          uxxpuxxm[r][s][m][n];
                     else
                        A[2+grid.dim+two2one[m][n]][i] += z[i][r]*z[i][s]*
                                                          uxxpuxxm[r][s][m][n];
            }

         if (globdebug)
         {
            tempval = A[0][i];
            tempval += A[1][i]*getu(x,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               tempval += A[2+m][i]*getDu(x,m,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  tempval += A[2+grid.dim+two2one[m][n]][i]*getD2u(x,m,n,thesign,grid);
            cout << "Taylor series diff = " 
                 << getu(index[i],0,0,0.0,-thesign,grid)-tempval << endl;

            tempval = 0.0;
            tempval += getu(x,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               tempval += z[i][m]*getDu(x,m,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  if (m == n)
                     tempval += 0.5*z[i][m]*z[i][n]*getD2u(x,m,n,thesign,grid);
                  else
                     tempval += z[i][m]*z[i][n]*getD2u(x,m,n,thesign,grid);
            cout << "   on other side 1 = " 
                 << getu(index[i],0,0,0.0,thesign,grid)-tempval << endl;

            tempval = 0.0;
            tempval += getu(x,-thesign,grid);
            for (m = 0; m < grid.dim; m++)
               tempval += z[i][m]*getDu(x,m,-thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  if (m == n)
                     tempval += 0.5*z[i][m]*z[i][n]*getD2u(x,m,n,-thesign,grid);
                  else
                     tempval += z[i][m]*z[i][n]*getD2u(x,m,n,-thesign,grid);
            cout << "   on other side 2 = " 
                 << getu(index[i],0,0,0.0,-thesign,grid)-tempval << endl;
         }

         double tempval2, tempval3, tempval4;
         if (globdebug)
         {
            tempval = 0.0;
            tempval2 = upum*getu(x,thesign,grid)+up0;
            cout << "      check: " << fabs(getu(x,-thesign,grid)-tempval2) << endl;
            temperr[0] = getu(x,-thesign,grid)-tempval2;
            tempval3 = temperr[0];
            tempval += tempval2;
            cout << "      more: " << getu(x,-thesign,grid)-tempval << endl;
            for (m = 0; m < grid.dim; m++)
            {
               tempval2 = uxpuxm[m][0]*getDu(x,0,thesign,grid)+
                          uxpuxm[m][1]*getDu(x,1,thesign,grid)+
                          uxpuxm[m][2]*getDu(x,2,thesign,grid)+uxp0[m];
               cout << "      check: " << fabs(getDu(x,m,-thesign,grid)-tempval2) << endl;
               temperr[1+m] = getDu(x,m,-thesign,grid)-tempval2;
               tempval3 += z[i][m]*temperr[1+m];
               tempval += z[i][m]*tempval2;
            }
            cout << "      more: " << getu(x,-thesign,grid)+
                                      z[i][0]*getDu(x,0,-thesign,grid)+
                                      z[i][1]*getDu(x,1,-thesign,grid)+
                                      z[i][2]*getDu(x,2,-thesign,grid)-tempval << endl;
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
               {
                  tempval2 = uxxpuxxm[m][n][0][0]*getD2u(x,0,0,thesign,grid)+
                             uxxpuxxm[m][n][0][1]*getD2u(x,0,1,thesign,grid)+
                             uxxpuxxm[m][n][0][2]*getD2u(x,0,2,thesign,grid)+
//                          uxxpuxxm[m][n][1][0]*getD2u(x,1,0,thesign,grid)+
                             uxxpuxxm[m][n][1][1]*getD2u(x,1,1,thesign,grid)+
                             uxxpuxxm[m][n][1][2]*getD2u(x,1,2,thesign,grid)+
//                          uxxpuxxm[m][n][2][0]*getD2u(x,2,0,thesign,grid)+
//                          uxxpuxxm[m][n][2][1]*getD2u(x,2,1,thesign,grid)+
                             uxxpuxxm[m][n][2][2]*getD2u(x,2,2,thesign,grid)+
                             uxxpuxm[m][n][0]*getDu(x,0,thesign,grid)+
                             uxxpuxm[m][n][1]*getDu(x,1,thesign,grid)+
                             uxxpuxm[m][n][2]*getDu(x,2,thesign,grid)+uxxp0[m][n];
                  cout << "      check: " << fabs(getD2u(x,m,n,-thesign,grid)-tempval2) 
                       << endl;
                  temperr[1+grid.dim+two2one[m][n]] = getD2u(x,m,n,-thesign,grid)-
                                                      tempval2;
                  if (m == n)
                     tempval3 += 0.5*z[i][m]*z[i][n]*temperr[1+grid.dim+two2one[m][n]];
                  else
                     tempval3 += z[i][m]*z[i][n]*temperr[1+grid.dim+two2one[m][n]];
                  if (m == n)
                     tempval += 0.5*z[i][m]*z[i][n]*tempval2;
                  else
                     tempval += z[i][m]*z[i][n]*tempval2;
               }
            cout << "      more: " << getu(x,-thesign,grid)+
                                      z[i][0]*getDu(x,0,-thesign,grid)+
                                      z[i][1]*getDu(x,1,-thesign,grid)+
                                      z[i][2]*getDu(x,2,-thesign,grid)+
                                      0.5*z[i][0]*z[i][0]*getD2u(x,0,0,-thesign,grid)+
                                      z[i][0]*z[i][1]*getD2u(x,0,1,-thesign,grid)+
                                      z[i][0]*z[i][2]*getD2u(x,0,2,-thesign,grid)+
                                      0.5*z[i][1]*z[i][1]*getD2u(x,1,1,-thesign,grid)+
                                      z[i][1]*z[i][2]*getD2u(x,1,2,-thesign,grid)+
                                      0.5*z[i][2]*z[i][2]*getD2u(x,2,2,-thesign,grid)-
                                      tempval << " " 
                                   << getu(x,-thesign,grid)+
                                      z[i][0]*getDu(x,0,-thesign,grid)+
                                      z[i][1]*getDu(x,1,-thesign,grid)+
                                      z[i][2]*getDu(x,2,-thesign,grid)+
                                      0.5*z[i][0]*z[i][0]*getD2u(x,0,0,-thesign,grid)+
                                      z[i][0]*z[i][1]*getD2u(x,0,1,-thesign,grid)+
                                      z[i][0]*z[i][2]*getD2u(x,0,2,-thesign,grid)+
                                      0.5*z[i][1]*z[i][1]*getD2u(x,1,1,-thesign,grid)+
                                      z[i][1]*z[i][2]*getD2u(x,1,2,-thesign,grid)+
                                      0.5*z[i][2]*z[i][2]*getD2u(x,2,2,-thesign,grid)-
                                      getu(index[i],0,0,0.0,-thesign,grid) << endl;
            cout << "   on other side again = " 
                 << getu(index[i],0,0,0.0,-thesign,grid)-tempval << endl;
            cout << "   on other side error = " << tempval3 << endl;
         }
      }
   }

   if (globdebug)
   {
      cout << "the matrix = " << endl;
      for (m = 0; m < M; m++)
      {
         if (yesC)
            for (n = 0; n <= N; n++)
               cout << A[m][n] << " ";
         else
            for (n = 0; n < N; n++)
               cout << A[m][n] << " ";
         cout << endl;
      }
   }

   free_matrix(z,N-1,grid.dim-1);
   delete [] signs;
}

void getiimCstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                      double up0, double upum, double *uxp0, double **uxpuxm, 
                      double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm, 
                      double ***S, PBData &pb, GridData &grid)
{

   double **z = matrix(N-1,grid.dim-1);
   int *signs = new int[N];
   int two2one[grid.dim][grid.dim];
   int i, r, s, m, n;

   cout << "main point " << x[0] << " " << x[1] << " " << x[2] << " has sign " 
        << thesign << endl;
   for (i = 0; i < N; i++)
   {
      sub2coord(z[i],index[i],grid);
      for (r = 0; r < grid.dim; r++)
         z[i][r] = z[i][r]-x[r];
      if (evalarray(S,index[i]) < 0.0)
         signs[i] = -1;
      else
         signs[i] = 1;
      cout << "point at " << index[i][0] << " " << index[i][1] << " " << index[i][2] 
           << " has sign " << signs[i] << " and dist " 
           << sqrt(getdotprod(z[i],z[i],grid.dim)) << endl;
   }
   for (n = 0; n < grid.dim; n++)
      two2one[n][n] = n;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (m = 0; m < grid.dim; m++)
      for (n = 0; n < m; n++)
         two2one[m][n] = two2one[n][m];

   double tempval, temperr[10];

   A[0][N] = 1.0;
   for (i = 1; i < M; i++)
      A[i][N] = 0.0;
   for (i = 0; i < N; i++)
   {
      if (signs[i] == thesign)
      {
         A[0][i] = 0.0;
         A[1][i] = 1.0;
         for (m = 0; m < grid.dim; m++)
            A[2+m][i] = z[i][m];
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               if (m == n)
                  A[2+grid.dim+two2one[m][n]][i] = 0.5*z[i][m]*z[i][n];
               else
                  A[2+grid.dim+two2one[m][n]][i] = z[i][m]*z[i][n];
         tempval = A[0][i];
         tempval += A[1][i]*getu(x,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            tempval += A[2+m][i]*getDu(x,m,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               tempval += A[2+grid.dim+two2one[m][n]][i]*getD2u(x,m,n,thesign,grid);
         cout << "Taylor series diff = " << getu(index[i],0,0,0.0,thesign,grid)-tempval 
              << endl;
      }
      else
      {
         A[0][i] = up0;
         for (r = 0; r < grid.dim; r++)
            A[0][i] += z[i][r]*uxp0[r];
         for (r = 0; r < grid.dim; r++)
            for (s = r; s < grid.dim; s++)
               if (r == s)
                  A[0][i] += 0.5*z[i][r]*z[i][s]*uxxp0[r][s];
               else
                  A[0][i] += z[i][r]*z[i][s]*uxxp0[r][s];

         A[1][i] = upum;

         for (m = 0; m < grid.dim; m++)
         {
            A[2+m][i] = 0.0;
            for (r = 0; r < grid.dim; r++)
               A[2+m][i] += z[i][r]*uxpuxm[r][m];
            for (r = 0; r < grid.dim; r++)
               for (s = r; s < grid.dim; s++)
                  if (r == s)
                     A[2+m][i] += 0.5*z[i][r]*z[i][s]*uxxpuxm[r][s][m];
                  else
                     A[2+m][i] += z[i][r]*z[i][s]*uxxpuxm[r][s][m];
         }

         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
            {
               A[2+grid.dim+two2one[m][n]][i] = 0.0;
               for (r = 0; r < grid.dim; r++)
                  for (s = r; s < grid.dim; s++)
                     if (r == s)
                        A[2+grid.dim+two2one[m][n]][i] += 0.5*z[i][r]*z[i][s]*
                                                          uxxpuxxm[r][s][m][n];
                     else
                        A[2+grid.dim+two2one[m][n]][i] += z[i][r]*z[i][s]*
                                                          uxxpuxxm[r][s][m][n];
            }

         tempval = A[0][i];
         tempval += A[1][i]*getu(x,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            tempval += A[2+m][i]*getDu(x,m,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               tempval += A[2+grid.dim+two2one[m][n]][i]*getD2u(x,m,n,thesign,grid);
         cout << "Taylor series diff = " << getu(index[i],0,0,0.0,-thesign,grid)-tempval 
              << endl;

         tempval = 0.0;
         tempval += getu(x,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            tempval += z[i][m]*getDu(x,m,thesign,grid);
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               if (m == n)
                  tempval += 0.5*z[i][m]*z[i][n]*getD2u(x,m,n,thesign,grid);
               else
                  tempval += z[i][m]*z[i][n]*getD2u(x,m,n,thesign,grid);
         cout << "   on other side 1 = " << getu(index[i],0,0,0.0,thesign,grid)-tempval 
              << endl;

         tempval = 0.0;
         tempval += getu(x,-thesign,grid);
         for (m = 0; m < grid.dim; m++)
            tempval += z[i][m]*getDu(x,m,-thesign,grid);
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               if (m == n)
                  tempval += 0.5*z[i][m]*z[i][n]*getD2u(x,m,n,-thesign,grid);
               else
                  tempval += z[i][m]*z[i][n]*getD2u(x,m,n,-thesign,grid);
         cout << "   on other side 2 = " << getu(index[i],0,0,0.0,-thesign,grid)-tempval 
              << endl;

         double tempval2, tempval3, tempval4;
         tempval = 0.0;
         tempval2 = upum*getu(x,thesign,grid)+up0;
         cout << "      check: " << fabs(getu(x,-thesign,grid)-tempval2) << endl;
         temperr[0] = getu(x,-thesign,grid)-tempval2;
         tempval3 = temperr[0];
         tempval += tempval2;
         cout << "      more: " << getu(x,-thesign,grid)-tempval << endl;
         for (m = 0; m < grid.dim; m++)
         {
            tempval2 = uxpuxm[m][0]*getDu(x,0,thesign,grid)+
                       uxpuxm[m][1]*getDu(x,1,thesign,grid)+
                       uxpuxm[m][2]*getDu(x,2,thesign,grid)+uxp0[m];
            cout << "      check: " << fabs(getDu(x,m,-thesign,grid)-tempval2) << endl;
            temperr[1+m] = getDu(x,m,-thesign,grid)-tempval2;
            tempval3 += z[i][m]*temperr[1+m];
            tempval += z[i][m]*tempval2;
         }
         cout << "      more: " << getu(x,-thesign,grid)+
                                   z[i][0]*getDu(x,0,-thesign,grid)+
                                   z[i][1]*getDu(x,1,-thesign,grid)+
                                   z[i][2]*getDu(x,2,-thesign,grid)-tempval << endl;
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
            {
               tempval2 = uxxpuxxm[m][n][0][0]*getD2u(x,0,0,thesign,grid)+
                          uxxpuxxm[m][n][0][1]*getD2u(x,0,1,thesign,grid)+
                          uxxpuxxm[m][n][0][2]*getD2u(x,0,2,thesign,grid)+
//                          uxxpuxxm[m][n][1][0]*getD2u(x,1,0,thesign,grid)+
                          uxxpuxxm[m][n][1][1]*getD2u(x,1,1,thesign,grid)+
                          uxxpuxxm[m][n][1][2]*getD2u(x,1,2,thesign,grid)+
//                          uxxpuxxm[m][n][2][0]*getD2u(x,2,0,thesign,grid)+
//                          uxxpuxxm[m][n][2][1]*getD2u(x,2,1,thesign,grid)+
                          uxxpuxxm[m][n][2][2]*getD2u(x,2,2,thesign,grid)+
                          uxxpuxm[m][n][0]*getDu(x,0,thesign,grid)+
                          uxxpuxm[m][n][1]*getDu(x,1,thesign,grid)+
                          uxxpuxm[m][n][2]*getDu(x,2,thesign,grid)+uxxp0[m][n];
               cout << "      check: " << fabs(getD2u(x,m,n,-thesign,grid)-tempval2) 
                    << endl;
               temperr[1+grid.dim+two2one[m][n]] = getD2u(x,m,n,-thesign,grid)-tempval2;
               if (m == n)
                  tempval3 += 0.5*z[i][m]*z[i][n]*temperr[1+grid.dim+two2one[m][n]];
               else
                  tempval3 += z[i][m]*z[i][n]*temperr[1+grid.dim+two2one[m][n]];
               if (m == n)
                  tempval += 0.5*z[i][m]*z[i][n]*tempval2;
               else
                  tempval += z[i][m]*z[i][n]*tempval2;
            }
         cout << "      more: " << getu(x,-thesign,grid)+
                                   z[i][0]*getDu(x,0,-thesign,grid)+
                                   z[i][1]*getDu(x,1,-thesign,grid)+
                                   z[i][2]*getDu(x,2,-thesign,grid)+
                                   0.5*z[i][0]*z[i][0]*getD2u(x,0,0,-thesign,grid)+
                                   z[i][0]*z[i][1]*getD2u(x,0,1,-thesign,grid)+
                                   z[i][0]*z[i][2]*getD2u(x,0,2,-thesign,grid)+
                                   0.5*z[i][1]*z[i][1]*getD2u(x,1,1,-thesign,grid)+
                                   z[i][1]*z[i][2]*getD2u(x,1,2,-thesign,grid)+
                                   0.5*z[i][2]*z[i][2]*getD2u(x,2,2,-thesign,grid)-
                                   tempval << " " 
                                << getu(x,-thesign,grid)+
                                   z[i][0]*getDu(x,0,-thesign,grid)+
                                   z[i][1]*getDu(x,1,-thesign,grid)+
                                   z[i][2]*getDu(x,2,-thesign,grid)+
                                   0.5*z[i][0]*z[i][0]*getD2u(x,0,0,-thesign,grid)+
                                   z[i][0]*z[i][1]*getD2u(x,0,1,-thesign,grid)+
                                   z[i][0]*z[i][2]*getD2u(x,0,2,-thesign,grid)+
                                   0.5*z[i][1]*z[i][1]*getD2u(x,1,1,-thesign,grid)+
                                   z[i][1]*z[i][2]*getD2u(x,1,2,-thesign,grid)+
                                   0.5*z[i][2]*z[i][2]*getD2u(x,2,2,-thesign,grid)-
                                   getu(index[i],0,0,0.0,-thesign,grid) << endl;
         cout << "   on other side again = " 
              << getu(index[i],0,0,0.0,-thesign,grid)-tempval << endl;
         cout << "   on other side error = " << tempval3 << endl;
      }
   }

   cout << "the matrix = " << endl;
   for (m = 0; m < M; m++)
   {
      for (n = 0; n <= N; n++)
         cout << A[m][n] << " ";
      cout << endl;
   }

   free_matrix(z,N-1,grid.dim-1);
   delete [] signs;
}

void getiimgridstencilmx(double **A, int M, int N, double *x, int *gindex, int thesign, 
                         int **index, char yesC, double up0, double upum, double *uxp0, 
                         double **uxpuxm, double **uxxp0, double ***uxxpuxm, 
                         double ****uxxpuxxm, double ***S, PBData &pb, GridData &grid)
{

   double **z = matrix(N-1,grid.dim-1), y[grid.dim];
   int *signs = new int[N];
   int two2one[grid.dim][grid.dim];
   int i, r, s, m, n;

   if (globdebug)
      cout << "main point " << x[0] << " " << x[1] << " " << x[2] << " has sign " 
           << thesign << endl;
   sub2coord(y,gindex,grid);
   for (n = 0; n < grid.dim; n++)
      two2one[n][n] = n;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (m = 0; m < grid.dim; m++)
      for (n = 0; n < m; n++)
         two2one[m][n] = two2one[n][m];

   double tempval, temperr[10];

   if (yesC)
   {
      A[0][N] = 1.0;
      for (i = 1; i < M; i++)
         A[i][N] = 0.0;
   }
   for (i = 0; i < N; i++)
   {
      sub2coord(z[i],index[i],grid);
      if (evalarray(S,index[i]) < 0.0)
         signs[i] = -1;
      else
         signs[i] = 1;
      if (globdebug)
         cout << "point at " << index[i][0] << " " << index[i][1] << " " << index[i][2] 
              << " has sign " << signs[i] << endl;
      if (signs[i] == thesign)
      {
         for (r = 0; r < grid.dim; r++)
            z[i][r] = z[i][r]-y[r];
         A[0][i] = 0.0;
         A[1][i] = 1.0;
         for (m = 0; m < grid.dim; m++)
            A[2+m][i] = z[i][m];
         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
               if (m == n)
                  A[2+grid.dim+two2one[m][n]][i] = 0.5*z[i][m]*z[i][n];
               else
                  A[2+grid.dim+two2one[m][n]][i] = z[i][m]*z[i][n];
         if (globdebug)
         {
            tempval = A[0][i];
            tempval += A[1][i]*getu(x,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               tempval += A[2+m][i]*getDu(x,m,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  tempval += A[2+grid.dim+two2one[m][n]][i]*getD2u(x,m,n,thesign,grid);
            cout << "Taylor series diff = " 
                 << getu(index[i],0,0,0.0,thesign,grid)-tempval << endl;
         }
      }
      else
      {
         for (r = 0; r < grid.dim; r++)
            z[i][r] = z[i][r]-x[r];
         A[0][i] = up0;
         for (r = 0; r < grid.dim; r++)
            A[0][i] += z[i][r]*uxp0[r];
         for (r = 0; r < grid.dim; r++)
            for (s = r; s < grid.dim; s++)
               if (r == s)
                  A[0][i] += 0.5*z[i][r]*z[i][s]*uxxp0[r][s];
               else
                  A[0][i] += z[i][r]*z[i][s]*uxxp0[r][s];

         A[1][i] = upum;

         for (m = 0; m < grid.dim; m++)
         {
            A[2+m][i] = 0.0;
            for (r = 0; r < grid.dim; r++)
               A[2+m][i] += z[i][r]*uxpuxm[r][m];
            for (r = 0; r < grid.dim; r++)
               for (s = r; s < grid.dim; s++)
                  if (r == s)
                     A[2+m][i] += 0.5*z[i][r]*z[i][s]*uxxpuxm[r][s][m];
                  else
                     A[2+m][i] += z[i][r]*z[i][s]*uxxpuxm[r][s][m];
         }

         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
            {
               A[2+grid.dim+two2one[m][n]][i] = 0.0;
               for (r = 0; r < grid.dim; r++)
                  for (s = r; s < grid.dim; s++)
                     if (r == s)
                        A[2+grid.dim+two2one[m][n]][i] += 0.5*z[i][r]*z[i][s]*
                                                          uxxpuxxm[r][s][m][n];
                     else
                        A[2+grid.dim+two2one[m][n]][i] += z[i][r]*z[i][s]*
                                                          uxxpuxxm[r][s][m][n];
            }

         for (r = 0; r < grid.dim; r++)
            z[i][r] = x[r]-y[r];

         for (m = 0; m < grid.dim; m++)
            for (n = m; n < grid.dim; n++)
            {
               if (m != n)
                  A[2+grid.dim+two2one[m][n]][i] += A[2+m][i]*z[i][n]+A[2+n][i]*z[i][m];
               else
                  A[2+grid.dim+two2one[m][n]][i] += A[2+m][i]*z[i][n];
               if (m != n)
                  A[2+grid.dim+two2one[m][n]][i] += A[1][i]*z[i][m]*z[i][n];
               else
                  A[2+grid.dim+two2one[m][n]][i] += 0.5*A[1][i]*z[i][m]*z[i][n];
            }
         for (m = 0; m < grid.dim; m++)
            A[2+m][i] += A[1][i]*z[i][m];

         if (globdebug)
         {
            tempval = A[0][i];
            tempval += A[1][i]*getu(x,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               tempval += A[2+m][i]*getDu(x,m,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  tempval += A[2+grid.dim+two2one[m][n]][i]*getD2u(x,m,n,thesign,grid);
            cout << "Taylor series diff = " 
                 << getu(index[i],0,0,0.0,-thesign,grid)-tempval << endl;

            tempval = 0.0;
            tempval += getu(x,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               tempval += z[i][m]*getDu(x,m,thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  if (m == n)
                     tempval += 0.5*z[i][m]*z[i][n]*getD2u(x,m,n,thesign,grid);
                  else
                     tempval += z[i][m]*z[i][n]*getD2u(x,m,n,thesign,grid);
            cout << "   on other side 1 = " 
                 << getu(index[i],0,0,0.0,thesign,grid)-tempval << endl;
   
            tempval = 0.0;
            tempval += getu(x,-thesign,grid);
            for (m = 0; m < grid.dim; m++)
               tempval += z[i][m]*getDu(x,m,-thesign,grid);
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
                  if (m == n)
                     tempval += 0.5*z[i][m]*z[i][n]*getD2u(x,m,n,-thesign,grid);
                  else
                     tempval += z[i][m]*z[i][n]*getD2u(x,m,n,-thesign,grid);
            cout << "   on other side 2 = " 
                 << getu(index[i],0,0,0.0,-thesign,grid)-tempval << endl;

            double tempval2, tempval3, tempval4;
            tempval = 0.0;
            tempval2 = upum*getu(x,thesign,grid)+up0;
            cout << "      check: " << fabs(getu(x,-thesign,grid)-tempval2) << endl;
            temperr[0] = getu(x,-thesign,grid)-tempval2;
            tempval3 = temperr[0];
            tempval += tempval2;
            cout << "      more: " << getu(x,-thesign,grid)-tempval << endl;
            for (m = 0; m < grid.dim; m++)
            {
               tempval2 = uxpuxm[m][0]*getDu(x,0,thesign,grid)+
                          uxpuxm[m][1]*getDu(x,1,thesign,grid)+
                          uxpuxm[m][2]*getDu(x,2,thesign,grid)+uxp0[m];
               cout << "      check: " << fabs(getDu(x,m,-thesign,grid)-tempval2) 
                    << endl;
               temperr[1+m] = getDu(x,m,-thesign,grid)-tempval2;
               tempval3 += z[i][m]*temperr[1+m];
               tempval += z[i][m]*tempval2;
            }
            cout << "      more: " << getu(x,-thesign,grid)+
                                      z[i][0]*getDu(x,0,-thesign,grid)+
                                      z[i][1]*getDu(x,1,-thesign,grid)+
                                      z[i][2]*getDu(x,2,-thesign,grid)-tempval << endl;
            for (m = 0; m < grid.dim; m++)
               for (n = m; n < grid.dim; n++)
               {
                  tempval2 = uxxpuxxm[m][n][0][0]*getD2u(x,0,0,thesign,grid)+
                             uxxpuxxm[m][n][0][1]*getD2u(x,0,1,thesign,grid)+
                             uxxpuxxm[m][n][0][2]*getD2u(x,0,2,thesign,grid)+
//                          uxxpuxxm[m][n][1][0]*getD2u(x,1,0,thesign,grid)+
                             uxxpuxxm[m][n][1][1]*getD2u(x,1,1,thesign,grid)+
                             uxxpuxxm[m][n][1][2]*getD2u(x,1,2,thesign,grid)+
//                          uxxpuxxm[m][n][2][0]*getD2u(x,2,0,thesign,grid)+
//                          uxxpuxxm[m][n][2][1]*getD2u(x,2,1,thesign,grid)+
                             uxxpuxxm[m][n][2][2]*getD2u(x,2,2,thesign,grid)+
                             uxxpuxm[m][n][0]*getDu(x,0,thesign,grid)+
                             uxxpuxm[m][n][1]*getDu(x,1,thesign,grid)+
                             uxxpuxm[m][n][2]*getDu(x,2,thesign,grid)+uxxp0[m][n];
                  cout << "      check: " << fabs(getD2u(x,m,n,-thesign,grid)-tempval2) 
                       << endl;
                  temperr[1+grid.dim+two2one[m][n]] = getD2u(x,m,n,-thesign,grid)-
                                                      tempval2;
                  if (m == n)
                     tempval3 += 0.5*z[i][m]*z[i][n]*temperr[1+grid.dim+two2one[m][n]];
                  else
                     tempval3 += z[i][m]*z[i][n]*temperr[1+grid.dim+two2one[m][n]];
                  if (m == n)
                     tempval += 0.5*z[i][m]*z[i][n]*tempval2;
                  else
                     tempval += z[i][m]*z[i][n]*tempval2;
               }
            cout << "      more: " << getu(x,-thesign,grid)+
                                      z[i][0]*getDu(x,0,-thesign,grid)+
                                      z[i][1]*getDu(x,1,-thesign,grid)+
                                      z[i][2]*getDu(x,2,-thesign,grid)+
                                      0.5*z[i][0]*z[i][0]*getD2u(x,0,0,-thesign,grid)+
                                      z[i][0]*z[i][1]*getD2u(x,0,1,-thesign,grid)+
                                      z[i][0]*z[i][2]*getD2u(x,0,2,-thesign,grid)+
                                      0.5*z[i][1]*z[i][1]*getD2u(x,1,1,-thesign,grid)+
                                      z[i][1]*z[i][2]*getD2u(x,1,2,-thesign,grid)+
                                      0.5*z[i][2]*z[i][2]*getD2u(x,2,2,-thesign,grid)-
                                      tempval << " " 
                                   << getu(x,-thesign,grid)+
                                      z[i][0]*getDu(x,0,-thesign,grid)+
                                      z[i][1]*getDu(x,1,-thesign,grid)+
                                      z[i][2]*getDu(x,2,-thesign,grid)+
                                      0.5*z[i][0]*z[i][0]*getD2u(x,0,0,-thesign,grid)+
                                      z[i][0]*z[i][1]*getD2u(x,0,1,-thesign,grid)+
                                      z[i][0]*z[i][2]*getD2u(x,0,2,-thesign,grid)+
                                      0.5*z[i][1]*z[i][1]*getD2u(x,1,1,-thesign,grid)+
                                      z[i][1]*z[i][2]*getD2u(x,1,2,-thesign,grid)+
                                      0.5*z[i][2]*z[i][2]*getD2u(x,2,2,-thesign,grid)-
                                      getu(index[i],0,0,0.0,-thesign,grid) << endl;
            cout << "   on other side again = " 
                 << getu(index[i],0,0,0.0,-thesign,grid)-tempval << endl;
            cout << "   on other side error = " << tempval3 << endl;
         }
      }
   }

   if (globdebug)
   {
      cout << "the matrix = " << endl;
      for (m = 0; m < M; m++)
      {
         if (yesC)
            for (n = 0; n <= N; n++)
               cout << A[m][n] << " ";
         else
            for (n = 0; n < N; n++)
               cout << A[m][n] << " ";
         cout << endl;
      }
   }

   free_matrix(z,N-1,grid.dim-1);
   delete [] signs;
}

void iim(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S,
         PBData &pb, GridData &grid)
{
   int M = 2+grid.dim+(grid.dim*(grid.dim+1))/2;
   int N = 1;
   int i, j, k, r, s, t, rstar, sstar;
   for (r = 0; r < grid.dim; r++)
      N *= 3;
   double **B = matrix(M-1,N-1), *d = new double[M], *c = new double[N];
   double **LU = matrix(M-1,N-1);
   int PLR[M], PLC[N];
   int tindex[grid.dim], **cindex = imatrix(N-1,grid.dim-1);
   double up0, upum, uxp0[grid.dim], **uxpuxm = matrix(grid.dim-1,grid.dim-1),
          **uxxp0 = matrix(grid.dim-1,grid.dim-1), 
          ***uxxpuxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1),
          ****uxxpuxxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1,grid.dim-1);
   double ehere, ethere, thesign, alpha, tempalpha; 
   double x[grid.dim], normal[grid.dim], tempnormal[grid.dim], junk[grid.dim];
   int two2one[grid.dim][grid.dim];

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
      for (s = -1; s <= 1; s += 2)
         if (gamma[r][(s+1)/2] == 1)
         {
            getinterfaceinfo(tempalpha,junk,tempnormal,S,index,r,s,grid);
            if (tempalpha < alpha)
            {
               alpha = tempalpha;
               rstar = r;
               sstar = s;
               for (t = 0; t < grid.dim; t++)
                  normal[t] = tempnormal[t];
               sub2coord(x,index,grid);
               x[rstar] += sstar*alpha*grid.dx[rstar];
            }
         }

   getiimjumps(up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm,index,rstar,sstar,alpha,
               thesign,normal,S,pb,grid);

   j = 0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = -1;
   while (tindex[0] <= 1)
   {
      for (k = 0; k < grid.dim; k++)
         cindex[j][k] = index[k]+tindex[k];
      j++;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > 1; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   getiimstencilmx(B,M,N,x,thesign,cindex,up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm, 
                   S,pb,grid);

   for (r = 0; r < grid.dim; r++)
      two2one[r][r] = r;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (s = 0; s < grid.dim; s++)
      for (r = 0; r < s; r++)
         two2one[s][r] = two2one[r][s];
   for (r = 0; r < M; r++)
      d[r] = 0.0;
   for (r = 0; r < grid.dim; r++)
      d[2+grid.dim+two2one[r][r]] = 1.0;

   gecp0(LU,PLR,PLC,B,M-1,N-1);
   forwardbacksub0(c,d,LU,PLR,PLC,M-1);

   for (r = 0; r < N; r++)
   {
      if (fabs(c[r]) > grid.tol)
         sparse2(index,cindex[r],A,c[r],grid);
   }
}

void iim(SparseElt2**** &A, double ***b, int *index, int add, double ***S, 
         char ***tube, PBData &pb, GridData &grid)
{
   int M = 2+grid.dim+(grid.dim*(grid.dim+1))/2;
   int N = 1;
   int i, j, k, r, s, t, rstar, sstar;
   for (r = 0; r < grid.dim; r++)
      N *= 3;
   double **B = matrix(N-1,N-1), *d = new double[N], *c = new double[N];
   double **LU = matrix(N-1,N-1);
   int PLR[N], PLC[N];
   int tindex[grid.dim], **cindex = imatrix(N,grid.dim-1);
   double up0, upum, uxp0[grid.dim], **uxpuxm = matrix(grid.dim-1,grid.dim-1),
          **uxxp0 = matrix(grid.dim-1,grid.dim-1), 
          ***uxxpuxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1),
          ****uxxpuxxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1,grid.dim-1);
   double ehere, ethere, thesign, alpha, tempalpha; 
   double x[grid.dim], normal[grid.dim], tempnormal[grid.dim], junk[grid.dim];
   int two2one[grid.dim][grid.dim];

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
      two2one[r][r] = r;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (s = 0; s < grid.dim; s++)
      for (r = 0; r < s; r++)
         two2one[s][r] = two2one[r][s];
  
   double tempuval[N];

   getnearestinterface(x,index,S,grid);
   getiimjumps(up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm,x,thesign,S,pb,grid);
   getnearestgrid(cindex,N,x,M+add,-1.0,tube,grid);
   cout << "using stencil: " << endl;
   for (r = 0; r < N; r++)
      cout << cindex[r][0] << " " << cindex[r][1] << " " << cindex[r][2] 
           << " with value " << evalarray(S,cindex[r]) << endl;
   cout << "number of points should be " << M << endl;
   getiimstencilmx(B,M,N,x,thesign,cindex,up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm, 
                   S,pb,grid);
//   getiimgridstencilmx(B,M,N,x,index,thesign,cindex,0,up0,upum,uxp0,uxpuxm,uxxp0,
//                       uxxpuxm,uxxpuxxm,S,pb,grid);
   for (r = 0; r < N; r++)
   {
      tempuval[r] = B[0][r];
      tempuval[r] += B[1][r]*getu(x,thesign,grid);
      for (s = 0; s < grid.dim; s++)
         tempuval[r] += B[s+2][r]*getDu(x,s,thesign,grid);
      for (s = 0; s < grid.dim; s++)
         for (t = s; t < grid.dim; t++)
            tempuval[r] += B[2+grid.dim+two2one[s][t]][r]*getD2u(x,s,t,thesign,grid);
      cout << "Taylor poly u is " << tempuval[r] << " and real u is ";
      if (evalarray(S,cindex[r]) < 0.0)
         cout << getu(cindex[r],0,0,0.0,-1,grid) << " and error is " 
              << fabs(getu(cindex[r],0,0,0.0,-1,grid)-tempuval[r]) << endl;
      else
         cout << getu(cindex[r],0,0,0.0,1,grid) << " and error is " 
              << fabs(getu(cindex[r],0,0,0.0,1,grid)-tempuval[r]) << endl;
   }

   for (r = 0; r < M; r++)
      d[r] = 0.0;
//   for (r = 0; r < grid.dim; r++)
//      d[2+grid.dim+two2one[r][r]] = 1.0;
   d[2+grid.dim+two2one[1][1]] = 1.0;
//   d[3] = 1.0;

   if (M != N)
   {
      for (r = 0; r < M; r++)
         for (s = 0; s < N; s++)
            LU[r][s] = B[r][s];
      for (r = 0; r < N; r++)
         for (s = r; s < N; s++)
         {
            B[r][s] = 0.0;
            for (t = 0; t < M; t++)
               B[r][s] += LU[t][r]*LU[t][s];
            B[s][r] = B[r][s];
         }
      for (r = 0; r < M; r++)
         c[r] = d[r];
      for (r = 0; r < N; r++)
      {
         d[r] = 0.0;
         for (t = 0; t < M; t++)
            d[r] += LU[t][r]*c[t];
      }
   }
   else
      cout << "these two should be the same: " << M << " and " << N << endl;

   gecp0(LU,PLR,PLC,B,N-1,N-1);
   forwardbacksub0(c,d,LU,PLR,PLC,N-1);

   double tempval;
   for (r = 0; r < M; r++)
   {
      tempval = d[r];
      for (s = 0; s < N; s++)
         tempval -= B[r][s]*c[s];
      cout << "res row " << r << " = " << tempval << " with d = " << d[r] << endl;
   }

   for (r = 0; r < N; r++)
   {
      if (fabs(c[r]) > grid.tol)
      {
         sparse2(index,cindex[r],A,c[r],grid);
         cout << c[r] << " " << cindex[r][0] << " " << cindex[r][1] << " " 
              << cindex[r][2] << endl;
      }
   }
//   setvalarray(b,index,evalarray(b,index)-c[N]);
//   cout << c[N] << endl;

   tempval = 0.0;
   for (r = 0; r < N; r++)
      if (evalarray(S,cindex[r]) < 0.0)
         tempval += c[r]*getu(cindex[r],0,0,0.0,-1,grid);
      else
         tempval += c[r]*getu(cindex[r],0,0,0.0,1,grid);
   cout << tempval << " " << getD2u(index,1,1,0,0,0.0,thesign,grid) << endl;
   cout << tempval << " " << getD2u(x,1,1,thesign,grid) << endl;
//   cout << tempval << " " << getDu(index,1,0,0,0.0,thesign,grid) << endl;
//   cout << tempval << " " << getDu(x,1,thesign,grid) << endl;
//   cout << tempval << " " << getu(index,0,0,0.0,thesign,grid) << endl;
//   cout << tempval << " " << getu(x,thesign,grid) << endl;

   double tempval2;
   tempval2 = 0.0;
   for (r = 0; r < N; r++)
      tempval2 += c[r]*tempuval[r];
   cout << tempval2 << " " << getD2u(x,1,1,thesign,grid) << endl;
//   cout << tempval2 << " " << getDu(x,1,thesign,grid) << endl;
//   cout << tempval2 << " " << getu(x,thesign,grid) << endl;

   double theerr[N];
   for (r = 0; r < N; r++)
   {
      tempval2 = c[r]*tempuval[r];
      if (evalarray(S,cindex[r]) < 0.0)
         tempval = c[r]*getu(cindex[r],0,0,0.0,-1,grid);
      else
         tempval = c[r]*getu(cindex[r],0,0,0.0,1,grid);
      theerr[r] = tempval-tempval2;
      cout << tempval << " " << tempval2 << " " << theerr[r] << endl;
   }

   free_matrix(B,N-1,N-1);
   free_matrix(LU,N-1,N-1);
   delete [] c;
   delete [] d;
}

void iim(SparseElt2**** &A, double ***b, int *index, char yesC, int add, double ***S, 
         char ***tube, PBData &pb, GridData &grid)
{
   int M = 2+grid.dim+(grid.dim*(grid.dim+1))/2;
   int N = 1;
   int i, j, k, r, s, t, rstar, sstar;
   for (r = 0; r < grid.dim; r++)
      N *= 3;
   int L;
   double **B, *d, *c, **LU, *w;
   int *PLR, *PLC;
   int tindex[grid.dim], **cindex; 
   double up0, upum, uxp0[grid.dim], **uxpuxm = matrix(grid.dim-1,grid.dim-1),
          **uxxp0 = matrix(grid.dim-1,grid.dim-1), 
          ***uxxpuxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1),
          ****uxxpuxxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1,grid.dim-1);
   double ehere, ethere, thesign, alpha, tempalpha; 
   double x[grid.dim], normal[grid.dim], tempnormal[grid.dim], junk[grid.dim];
   int two2one[grid.dim][grid.dim];

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
      two2one[r][r] = r;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (s = 0; s < grid.dim; s++)
      for (r = 0; r < s; r++)
         two2one[s][r] = two2one[r][s];
  
   double tempuval[N];

   cout << "in here" << endl;
   getnearestinterface(x,index,S,grid);
//   x[0] = -0.4622331906382491; 
//   x[1] = -0.1849393522817337; 
//   x[2] = -0.0462368693065115;
//   thesign = 1.0;
//   x[0] = 0.2001;
//   x[1] = 0.001;
//   x[2] = 0.001;
//   thesign = -1.0;
//   printf("%4.16f %4.16f %4.16f\n",x[0],x[1],x[2]);
//   cout << thesign << endl;
//   index[0] = (x[0]+1.0)/grid.dx[0];
//   index[1] = (x[1]+1.0)/grid.dx[1];
//   index[2] = (x[2]+1.0)/grid.dx[2];
   double y[grid.dim];
   sub2coord(y,index,grid);
   if (globdebug)
      cout << "dist from gridpt to interfacept = " << getdist(x,y,grid.dim) << endl;
   cout << "at here" << endl;
   getiimjumps(up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm,x,thesign,S,pb,grid);
   if (yesC)
   {
      cindex = imatrix(M-1+add-1,grid.dim-1);
      getnearestgrid(cindex,N,x,M-1+add,-1.0,tube,grid);
      L = N+1;
   }
   else
   {
      cindex = imatrix(M+add-1,grid.dim-1);
      getnearestgrid(cindex,N,x,M+add,-1.0,tube,grid);
      L = N;
   }

/*
   cindex = imatrix(N-1,grid.dim-1);
   L = 0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = -1;
   while (tindex[0] <= 1)
   {
      for (k = 0; k < grid.dim; k++)
         cindex[L][k] = index[k]+tindex[k];
      if (globdebug)
         cout << L << "  " << cindex[L][0] << " " << cindex[L][1] << " " << cindex[L][2] 
                   << "  " << tindex[0] << " " << tindex[1] << " " << tindex[2] << endl;
      L++;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > 1; i--)
      {
         tindex[i] = -1;
         (tindex[i-1])++;
      }
   }
*/
   if (globdebug)
      cout << "number of points selected = " << L << endl;

   B = matrix(L-1,L-1); 
   d = new double[L]; 
   c = new double[L];
   LU = matrix(L-1,L-1);
   PLR = new int[L]; 
   PLC = new int[L];
   w = new double[L];

   if (globdebug)
   {
      cout << "using stencil: " << endl;
      for (r = 0; r < N; r++)
         cout << cindex[r][0] << " " << cindex[r][1] << " " << cindex[r][2] 
              << " with value " << evalarray(S,cindex[r]) << endl;
      cout << "number of points should be " << M << endl;
   }
//   getiimstencilmx(B,M,N,x,thesign,cindex,up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm, 
//                   S,pb,grid);
   cout << "down here" << endl;
   getiimstencilmx(B,M,N,x,thesign,cindex,yesC,up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,
                   uxxpuxxm,S,pb,grid);
//   getiimgridstencilmx(B,M,N,x,index,thesign,cindex,yesC,up0,upum,uxp0,uxpuxm,uxxp0,
//                       uxxpuxm,uxxpuxxm,S,pb,grid);
   if (globdebug)
      for (r = 0; r < N; r++)
      {
         tempuval[r] = B[0][r];
         tempuval[r] += B[1][r]*getu(x,thesign,grid);
         for (s = 0; s < grid.dim; s++)
            tempuval[r] += B[s+2][r]*getDu(x,s,thesign,grid);
         for (s = 0; s < grid.dim; s++)
            for (t = s; t < grid.dim; t++)
               tempuval[r] += B[2+grid.dim+two2one[s][t]][r]*getD2u(x,s,t,thesign,grid);
         cout << "Taylor poly u is " << tempuval[r] << " and real u is ";
         if (evalarray(S,cindex[r]) < 0.0)
            cout << getu(cindex[r],0,0,0.0,-1,grid) << " and error is " 
                 << fabs(getu(cindex[r],0,0,0.0,-1,grid)-tempuval[r]) << endl;
         else
            cout << getu(cindex[r],0,0,0.0,1,grid) << " and error is " 
                 << fabs(getu(cindex[r],0,0,0.0,1,grid)-tempuval[r]) << endl;
/*
         if (r == 1)
         {
            cout << getu(cindex[r],0,0,0.0,1,grid) << " " 
                 << getu(cindex[r],0,0,0.0,1,grid) << endl;
            cout << getu(x,1,grid)+grid.dx[0]*getDu(x,0,1,grid)+
                    0.5*grid.dx[0]*grid.dx[0]*getD2u(x,0,0,1,grid) <<  " "
                 << getu(x,-1,grid)+grid.dx[0]*getDu(x,0,-1,grid)+
                    0.5*grid.dx[0]*grid.dx[0]*getD2u(x,0,0,-1,grid) <<  endl;
         }
*/
      }

   for (r = 0; r < M; r++)
      d[r] = 0.0;
   for (r = 0; r < grid.dim; r++)
      if (thesign < 0.0)
         d[2+grid.dim+two2one[r][r]] = -pb.epsilonm;
      else
         d[2+grid.dim+two2one[r][r]] = -pb.epsilonp;
//      d[2+grid.dim+two2one[r][r]] = 1.0;
//   d[2+grid.dim+two2one[0][0]] = 1.0;
//   d[1] = 1.0;

   if (M != L)
   {
//      svdcmp(B,M,L,w,LU);
//      svbksb(B,w,LU,M,L,d,c,1.0e-14);
      double **AA = matrix(M-1,L-1), dd[M];
      for (r = 0; r < M; r++)
      {
         for (s = 0; s < L; s++)
            AA[r][s] = B[r][s];
         dd[r] = d[r];
      }
      cout << "penalty here " << N << endl;
      double penalty[N], maxpenalty = 1.0e8, minpenalty = 1.0e-8;
      for (s = 0; s < N; s++)
      {
         sub2coord(y,cindex[s],grid);
         penalty[s] = getdist(x,y,grid.dim);
         if (penalty[s] > maxpenalty)
            penalty[s] = maxpenalty;
         if (penalty[s] < minpenalty)
            penalty[s] = minpenalty;
//         penalty[s] = 1.0;
         for (r = 0; r < M; r++)
            AA[r][s] /= penalty[s];
      }
//      t = 0;
//      for (s = 0; s < N && t != grid.dim; s++)
//         for (t = 0; t < grid.dim && cindex[s][t] == index[t]; t++);
//      s--;
//      for (r = 0; r < M; r++)
//         AA[r][s] /= penalty[s];
      
      svdcmp(AA,M,L,w,LU);
      svbksb(AA,w,LU,M,L,dd,c,1.0e-14);
      for (s = 0; s < N; s++)
         c[s] /= penalty[s];
      free_matrix(AA,M-1,L-1);
      if (globdebug)
      {
         cout << "soln = ";
         for (r = 0; r < L; r++)
            cout << c[r] << " ";
         cout << endl;
/*
         for (r = 0; r < L; r++)
         {
            c[r] = 0.0;
            for (s = 0; s < grid.dim && cindex[r][s] == index[s]; s++);
            if (s == grid.dim)
               c[r] = 2.0*grid.dim/(grid.dx[0]*grid.dx[0]);
            for (s = 0; s < grid.dim; s++)
               tindex[s] = index[s];
            for (i = 0; i < grid.dim; i++)
               for (j = -1; j <= 1; j += 2)
               {
                  tindex[i] = index[i]+j;
                  for (s = 0; s < grid.dim && cindex[r][s] == tindex[s]; s++);
                  if (s == grid.dim)
                     c[r] = -1.0/(grid.dx[0]*grid.dx[0]);
                  tindex[i] = index[i];
               }
         }
*/
      }   
   }
   else
   {
      if (globdebug)
         cout << "these two should be the same: " << M << " and " << L << endl;
      gecp0(LU,PLR,PLC,B,L-1,L-1);
      forwardbacksub0(c,d,LU,PLR,PLC,L-1);
   }

   double tempval;
   if (globdebug)
   {
      cout << "matrix is " << M << "x" << L << endl;
      for (r = 0; r < M; r++)
      {
         tempval = d[r];
         for (s = 0; s < L; s++)
            tempval -= B[r][s]*c[s];
         cout << "res row " << r << " = " << tempval << " with d = " << d[r] << endl;
      }
   }

   cout << "sparse here" << endl;
   for (r = 0; r < N; r++)
   {
      if (fabs(c[r]) > grid.tol)
      {
         sparse2(index,cindex[r],A,c[r],grid);
         if (globdebug)
            cout << c[r] << " " << cindex[r][0] << " " << cindex[r][1] << " " 
                 << cindex[r][2] << endl;
      }
   }
   if (yesC)
   {
      setvalarray(b,index,evalarray(b,index)-c[N]);
      if (globdebug)
         cout << c[N] << endl;
   }

   if (globdebug)
   {
      if (yesC)
         tempval = c[N];
      else
         tempval = 0.0;
      for (r = 0; r < N; r++)
         if (evalarray(S,cindex[r]) < 0.0)
            tempval += c[r]*getu(cindex[r],0,0,0.0,-1,grid);
         else
            tempval += c[r]*getu(cindex[r],0,0,0.0,1,grid);
      if (thesign < 0.0) 
      {
         cout << tempval << " " << -pb.epsilonm*
                                   (getD2u(index,0,0,0,0,0.0,thesign,grid)+ 
                                    getD2u(index,1,1,0,0,0.0,thesign,grid)+ 
                                    getD2u(index,2,2,0,0,0.0,thesign,grid)) << endl;
         cout << tempval << " " << -pb.epsilonm*
                                   (getD2u(x,0,0,thesign,grid)+ 
                                    getD2u(x,1,1,thesign,grid)+ 
                                    getD2u(x,2,2,thesign,grid)) << endl;
         cout << -pb.epsilonm*(getD2u(x,0,0,thesign,grid)+
                               getD2u(x,1,1,thesign,grid)+
                               getD2u(x,2,2,thesign,grid))-tempval << endl;
         cout << "   " << getD2u(x,0,0,thesign,grid) << " "
                       << getD2u(x,1,1,thesign,grid) << " "
                       << getD2u(x,2,2,thesign,grid) << " " << pb.epsilonm << endl;
      }
      else
      {
         cout << tempval << " " << -pb.epsilonp*
                                   (getD2u(index,0,0,0,0,0.0,thesign,grid)+ 
                                    getD2u(index,1,1,0,0,0.0,thesign,grid)+ 
                                    getD2u(index,2,2,0,0,0.0,thesign,grid)) << endl;
         cout << tempval << " " << -pb.epsilonp*
                                   (getD2u(x,0,0,thesign,grid)+ 
                                    getD2u(x,1,1,thesign,grid)+ 
                                    getD2u(x,2,2,thesign,grid)) << endl;
      }
//      cout << tempval << " " << (getD2u(index,0,0,0,0,0.0,thesign,grid)+ 
//                                 getD2u(index,1,1,0,0,0.0,thesign,grid)+ 
//                                 getD2u(index,2,2,0,0,0.0,thesign,grid)) << endl;
//      cout << tempval << " " << (getD2u(x,0,0,thesign,grid)+ 
//                                 getD2u(x,1,1,thesign,grid)+ 
//                                 getD2u(x,2,2,thesign,grid)) << endl;
//      cout << tempval << " " << getD2u(index,0,0,0,0,0.0,thesign,grid) << endl;
//      cout << tempval << " " << getD2u(x,0,0,thesign,grid) << endl;
//      cout << tempval << " " << getu(x,thesign,grid) << endl;
   }
//   cout << tempval << " " << getDu(index,1,0,0,0.0,thesign,grid) << endl;
//   cout << tempval << " " << getDu(x,1,thesign,grid) << endl;
//   cout << tempval << " " << getu(index,0,0,0.0,thesign,grid) << endl;

   double tempval2;
   if (globdebug)
   {
      if (yesC)
         tempval2 = c[N];
      else
         tempval2 = 0.0;
      for (r = 0; r < N; r++)
         tempval2 += c[r]*tempuval[r];
      if (thesign < 0.0)
         cout << tempval2 << " " << -pb.epsilonm*
                                    (getD2u(index,0,0,0,0,0.0,thesign,grid)+ 
                                     getD2u(index,1,1,0,0,0.0,thesign,grid)+ 
                                     getD2u(index,2,2,0,0,0.0,thesign,grid)) << endl;
      else
         cout << tempval2 << " " << -pb.epsilonp*
                                    (getD2u(index,0,0,0,0,0.0,thesign,grid)+ 
                                     getD2u(index,1,1,0,0,0.0,thesign,grid)+ 
                                     getD2u(index,2,2,0,0,0.0,thesign,grid)) << endl;
//      cout << tempval2 << " " << (getD2u(index,0,0,0,0,0.0,thesign,grid)+ 
//                                  getD2u(index,1,1,0,0,0.0,thesign,grid)+ 
//                                  getD2u(index,2,2,0,0,0.0,thesign,grid)) << endl;
//      cout << tempval2 << " " << getu(x,thesign,grid) << endl;
   }
//   cout << tempval2 << " " << getDu(x,1,thesign,grid) << endl;

   double theerr[L];
   if (globdebug)
   {
      if (yesC)
      {
         tempval2 = c[N];
         tempval = c[N];
         theerr[N] = tempval-tempval2;
         cout << tempval << " " << tempval2 << " " << theerr[N] << endl;
      }
      for (r = 0; r < N; r++)
      {
         tempval2 = c[r]*tempuval[r];
         if (evalarray(S,cindex[r]) < 0.0)
            tempval = c[r]*getu(cindex[r],0,0,0.0,-1,grid);
         else
            tempval = c[r]*getu(cindex[r],0,0,0.0,1,grid);
         theerr[r] = tempval-tempval2;
         cout << tempval << " " << tempval2 << " " << theerr[r] << endl;
      }
   }

   if (globdebug)
   {
      cout << "all done" << endl;
      exit(1);
   }

   cout << "done here" << endl;
   free_matrix(B,L-1,L-1);
   free_matrix(LU,L-1,L-1);
   free_matrix(cindex,N-1,grid.dim-1);
   delete [] c;
   delete [] d;
   delete [] PLR;
   delete [] PLC;
   delete [] w;
   cout << "finished here" << endl;
}

void iimghost(SparseElt2**** &A, double ***b, int *index, char yesC, int add, 
              double ***S, char ***tube, PBData &pb, GridData &grid)
{
   int M = 2+grid.dim+(grid.dim*(grid.dim+1))/2;
   int N = 1;
   int i, j, k, r, s, t, m, n, rstar, sstar;
   for (r = 0; r < grid.dim; r++)
      N *= 3;
   int L;
   double **B, *d, *c, **LU, *w;
   int *PLR, *PLC;
   int tindex[grid.dim], **cindex; 
   double up0, upum, uxp0[grid.dim], **uxpuxm = matrix(grid.dim-1,grid.dim-1),
          **uxxp0 = matrix(grid.dim-1,grid.dim-1), 
          ***uxxpuxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1),
          ****uxxpuxxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1,grid.dim-1);
   double ehere, ethere, thesign, alpha, tempalpha; 
   double x[grid.dim], normal[grid.dim], tempnormal[grid.dim], junk[grid.dim];
   int two2one[grid.dim][grid.dim];
   char yesghost = 0;

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
      two2one[r][r] = r;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (s = 0; s < grid.dim; s++)
      for (r = 0; r < s; r++)
         two2one[s][r] = two2one[r][s];
  
   for (r = 0; r < grid.dim; r++)
      tindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
   {
      double mainalpha = 0.0;
      for (s = -1; s <= 1; s += 2)
      {
         tindex[r] = index[r]+s;
         if (evalarray(S,tindex)*evalarray(S,index) < 0.0)
            getinterfaceinfo(alpha,junk,junk,S,index,tindex,grid);
         else
            alpha = 1.0;
//         cout << "r, s, alpha = " << r << " " << s << " " << alpha << endl;
         mainalpha += alpha;
         tindex[r] = index[r];
      }
//      cout << "r, mainalpha = " << r << " " << mainalpha << endl;
      for (s = -1; s <= 1; s += 2)
      {
         tindex[r] = index[r]+s;
         if (evalarray(S,tindex)*evalarray(S,index) < 0.0)
         {
            getinterfaceinfo(alpha,junk,junk,S,index,tindex,grid);
            if (yesghost)
               sparse2(index,index,A,ehere/(alpha*grid.dx[r]*grid.dx[r]),grid);
            else
               sparse2(index,index,A,ehere/(alpha*grid.dx[r]*0.5*mainalpha*grid.dx[r]),
                       grid);
            sub2coord(x,index,grid);
            x[r] += alpha*s*grid.dx[r];
            getiimjumps(up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm,x,thesign,S,pb,grid);
            if (yesC)
            {
               cindex = imatrix(M-1+add-1,grid.dim-1);
               getnearestgrid(cindex,N,x,M-1+add,-1.0,tube,grid);
               L = N+1;
            }
            else
            {
               cindex = imatrix(M+add-1,grid.dim-1);
               getnearestgrid(cindex,N,x,M+add,-1.0,tube,grid);
               L = N;
            }

            B = matrix(L-1,L-1); 
            d = new double[L]; 
            c = new double[L];
            LU = matrix(L-1,L-1);
            PLR = new int[L]; 
            PLC = new int[L];
            w = new double[L];

            getiimstencilmx(B,M,N,x,thesign,cindex,yesC,up0,upum,uxp0,uxpuxm,uxxp0,
                            uxxpuxm,uxxpuxxm,S,pb,grid);
//            getiimgridstencilmx(B,M,N,x,index,thesign,cindex,yesC,up0,upum,uxp0,
//                                uxpuxm,uxxp0,uxxpuxm,uxxpuxxm,S,pb,grid);
            for (m = 0; m < M; m++)
               d[m] = 0.0;
            d[1] = 1.0;
            double **AA = matrix(M-1,L-1), dd[M];
            for (m = 0; m < M; m++)
            {
               for (n = 0; n < L; n++)
                  AA[m][n] = B[m][n];
               dd[m] = d[m];
            }
//            double smallnum = 0.000001;
//            t = 0;
//            for (n = 0; n < N && t != grid.dim; n++)
//               for (t = 0; t < grid.dim && cindex[n][t] == index[t]; t++);
//            n--;
//            for (m = 0; m < M; m++)
//               AA[m][n] /= smallnum;
            svdcmp(AA,M,L,w,LU);
            svbksb(AA,w,LU,M,L,dd,c,1.0e-14);
//            c[n] /= smallnum;
            free_matrix(AA,M-1,L-1);
            double tempval = 0.0;
            for (n = 0; n < N; n++)
            {
               double y[grid.dim];
               sub2coord(y,cindex[n],grid);
               if (evalarray(S,cindex[n]) < 0.0)
                  tempval += c[n]*getu(y,-1.0,grid);
               else
                  tempval += c[n]*getu(y,1.0,grid);
//               cout << cindex[n][0] << " " << cindex[n][1] << " " 
//                    << cindex[n][2] << "  " << tempval << endl;
            }
//            cout << "compare " << tempval << " " << getu(x,thesign,grid) << endl;
            for (n = 0; n < N; n++)
               if (fabs(c[n]) > grid.tol)
               {
                  if (yesghost)
                     sparse2(index,cindex[n],A,-ehere*c[n]/(alpha*grid.dx[r]*grid.dx[r]),
                             grid);
                  else
                     sparse2(index,cindex[n],A,-ehere*c[n]/(alpha*grid.dx[r]*0.5*
                                                            mainalpha*grid.dx[r]),grid);
//                  sparse2(index,cindex[n],A,-ehere*c[n],grid);
               }
            if (yesC)
               setvalarray(b,index,evalarray(b,index)+
                                   ehere*c[N]/(alpha*grid.dx[r]*0.5*mainalpha*
                                               grid.dx[r]));

            free_matrix(B,L-1,L-1);
            free_matrix(LU,L-1,L-1);
            free_matrix(cindex,N-1,grid.dim-1);
            delete [] c;
            delete [] d;
            delete [] PLR;
            delete [] PLC;
            delete [] w;
         }
         else
         {
            if (yesghost)
            {
               sparse2(index,tindex,A,-ehere/(grid.dx[r]*grid.dx[r]),grid);
               sparse2(index,index,A,ehere/(grid.dx[r]*grid.dx[r]),grid);
            }
            else
            {
               sparse2(index,tindex,A,-ehere/(grid.dx[r]*0.5*mainalpha*grid.dx[r]),grid);
               sparse2(index,index,A,ehere/(grid.dx[r]*0.5*mainalpha*grid.dx[r]),grid);
            }
         }
         tindex[r] = index[r];
      }
   }
/*
   SparseElt2 *current;
   double tempval = 0.0, y[grid.dim];
   for (current = evalarray(A,index); current != NULL; 
        current = (*current).next)
   {
      sub2coord(y,(*current).cindex,grid);
      if (evalarray(S,(*current).cindex) < 0.0)
         tempval += (*current).val*getu(y,-1.0,grid);
      else
         tempval += (*current).val*getu(y,1.0,grid);
//      cout << (*current).cindex[0] << " " << (*current).cindex[1] << " " 
//           << (*current).cindex[2] << "  " << (*current).val << endl;
   }
   sub2coord(x,index,grid);
   cout << "final val = " << tempval << " " << -ehere*(getD2u(x,0,0,thesign,grid)+
                                                       getD2u(x,1,1,thesign,grid)+
                                                       getD2u(x,2,2,thesign,grid)) 
        << " at " << index[0] << " " << index[1] << " " << index[2] << endl;
   getchar();
*/
}

void iimC(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S,
          PBData &pb, GridData &grid)
{
   int M = 2+grid.dim+(grid.dim*(grid.dim+1))/2;
   int N = 1;
   int i, j, k, r, s, t, rstar, sstar;
   for (r = 0; r < grid.dim; r++)
      N *= 3;
   double **B = matrix(M-1,N), *d = new double[M], *c = new double[N+1];
   double **LU = matrix(M-1,N);
   int PLR[M], PLC[N+1];
   int tindex[grid.dim], **cindex = imatrix(N,grid.dim-1);
   double up0, upum, uxp0[grid.dim], **uxpuxm = matrix(grid.dim-1,grid.dim-1),
          **uxxp0 = matrix(grid.dim-1,grid.dim-1), 
          ***uxxpuxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1),
          ****uxxpuxxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1,grid.dim-1);
   double ehere, ethere, thesign, alpha, tempalpha; 
   double x[grid.dim], normal[grid.dim], tempnormal[grid.dim], junk[grid.dim];
   int two2one[grid.dim][grid.dim];

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
      for (s = -1; s <= 1; s += 2)
         if (gamma[r][(s+1)/2] == 1)
         {
            getinterfaceinfo(tempalpha,junk,tempnormal,S,index,r,s,grid);
            if (tempalpha < alpha)
            {
               alpha = tempalpha;
               rstar = r;
               sstar = s;
               for (t = 0; t < grid.dim; t++)
                  normal[t] = tempnormal[t];
               sub2coord(x,index,grid);
               x[rstar] += sstar*alpha*grid.dx[rstar];
            }
         }

   getiimjumps(up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm,index,rstar,sstar,alpha,
               thesign,normal,S,pb,grid);

   j = 0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = -1;
   while (tindex[0] <= 1)
   {
      for (k = 0; k < grid.dim; k++)
         cindex[j][k] = index[k]+tindex[k];
      j++;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > 1; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   getiimCstencilmx(B,M,N,x,thesign,cindex,up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm, 
                    S,pb,grid);

   for (r = 0; r < grid.dim; r++)
      two2one[r][r] = r;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (s = 0; s < grid.dim; s++)
      for (r = 0; r < s; r++)
         two2one[s][r] = two2one[r][s];
   for (r = 0; r < M; r++)
      d[r] = 0.0;
   for (r = 0; r < grid.dim; r++)
      d[2+grid.dim+two2one[r][r]] = 1.0;

   gecp0(LU,PLR,PLC,B,M-1,N);
   forwardbacksub0(c,d,LU,PLR,PLC,M-1);

   for (r = 0; r < N; r++)
   {
      if (fabs(c[r]) > grid.tol)
         sparse2(index,cindex[r],A,c[r],grid);
   }
   setvalarray(b,index,evalarray(b,index)-c[N]);
}

void iimC(SparseElt2**** &A, double ***b, int *index, int add, double ***S, 
          char ***tube, PBData &pb, GridData &grid)
{
   int M = 2+grid.dim+(grid.dim*(grid.dim+1))/2;
   int N = 1;
   int i, j, k, r, s, t, rstar, sstar;
   for (r = 0; r < grid.dim; r++)
      N *= 3;
   double **B = matrix(N,N), *d = new double[N+1], *c = new double[N+1];
   double **LU = matrix(N,N);
   int PLR[N+1], PLC[N+1];
   int tindex[grid.dim], **cindex = imatrix(N,grid.dim-1);
   double up0, upum, uxp0[grid.dim], **uxpuxm = matrix(grid.dim-1,grid.dim-1),
          **uxxp0 = matrix(grid.dim-1,grid.dim-1), 
          ***uxxpuxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1),
          ****uxxpuxxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1,grid.dim-1);
   double ehere, ethere, thesign, alpha, tempalpha; 
   double x[grid.dim], normal[grid.dim], tempnormal[grid.dim], junk[grid.dim];
   int two2one[grid.dim][grid.dim];

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
      two2one[r][r] = r;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (s = 0; s < grid.dim; s++)
      for (r = 0; r < s; r++)
         two2one[s][r] = two2one[r][s];
  
   double tempuval[N];

   getnearestinterface(x,index,S,grid);
   getiimjumps(up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm,x,thesign,S,pb,grid);
   getnearestgrid(cindex,N,x,M-1+add,-1.0,tube,grid);
   cout << "using stencil: " << endl;
   for (r = 0; r < N; r++)
      cout << cindex[r][0] << " " << cindex[r][1] << " " << cindex[r][2] 
           << " with value " << evalarray(S,cindex[r]) << endl;
   cout << "number of points should be 1 less than " << M << endl;
   getiimCstencilmx(B,M,N,x,thesign,cindex,up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm, 
                    S,pb,grid);
//   getiimgridstencilmx(B,M,N,x,index,thesign,cindex,1,up0,upum,uxp0,uxpuxm,uxxp0,
//                       uxxpuxm,uxxpuxxm,S,pb,grid);
   for (r = 0; r < N; r++)
   {
      tempuval[r] = B[0][r];
      tempuval[r] += B[1][r]*getu(x,thesign,grid);
      for (s = 0; s < grid.dim; s++)
         tempuval[r] += B[s+2][r]*getDu(x,s,thesign,grid);
      for (s = 0; s < grid.dim; s++)
         for (t = s; t < grid.dim; t++)
            tempuval[r] += B[2+grid.dim+two2one[s][t]][r]*getD2u(x,s,t,thesign,grid);
      cout << "Taylor poly u is " << tempuval[r] << " and real u is ";
      if (evalarray(S,cindex[r]) < 0.0)
         cout << getu(cindex[r],0,0,0.0,-1,grid) << " and error is " 
              << fabs(getu(cindex[r],0,0,0.0,-1,grid)-tempuval[r]) << endl;
      else
         cout << getu(cindex[r],0,0,0.0,1,grid) << " and error is " 
              << fabs(getu(cindex[r],0,0,0.0,1,grid)-tempuval[r]) << endl;
   }

   for (r = 0; r < M; r++)
      d[r] = 0.0;
//   for (r = 0; r < grid.dim; r++)
//      d[2+grid.dim+two2one[r][r]] = 1.0;
   d[2+grid.dim+two2one[1][1]] = 1.0;
//   d[3] = 1.0;

   if (M-1 != N)
   {
      for (r = 0; r < M; r++)
         for (s = 0; s <= N; s++)
            LU[r][s] = B[r][s];
      for (r = 0; r <= N; r++)
         for (s = r; s <= N; s++)
         {
            B[r][s] = 0.0;
            for (t = 0; t < M; t++)
               B[r][s] += LU[t][r]*LU[t][s];
            B[s][r] = B[r][s];
         }
      for (r = 0; r < M; r++)
         c[r] = d[r];
      for (r = 0; r <= N; r++)
      {
         d[r] = 0.0;
         for (t = 0; t < M; t++)
            d[r] += LU[t][r]*c[t];
      }
   }
   else
      cout << "these two should be the same: " << M-1 << " and " << N << endl;

   gecp0(LU,PLR,PLC,B,N,N);
   forwardbacksub0(c,d,LU,PLR,PLC,N);

   double tempval;
   for (r = 0; r < M; r++)
   {
      tempval = d[r];
      for (s = 0; s <= N; s++)
         tempval -= B[r][s]*c[s];
      cout << "res row " << r << " = " << tempval << " with d = " << d[r] << endl;
   }

   for (r = 0; r < N; r++)
   {
      if (fabs(c[r]) > grid.tol)
      {
         sparse2(index,cindex[r],A,c[r],grid);
         cout << c[r] << " " << cindex[r][0] << " " << cindex[r][1] << " " 
              << cindex[r][2] << endl;
      }
   }
   setvalarray(b,index,evalarray(b,index)-c[N]);
   cout << c[N] << endl;

   tempval = c[N];
   for (r = 0; r < N; r++)
      if (evalarray(S,cindex[r]) < 0.0)
         tempval += c[r]*getu(cindex[r],0,0,0.0,-1,grid);
      else
         tempval += c[r]*getu(cindex[r],0,0,0.0,1,grid);
   cout << tempval << " " << getD2u(index,1,1,0,0,0.0,thesign,grid) << endl;
   cout << tempval << " " << getD2u(x,1,1,thesign,grid) << endl;
//   cout << tempval << " " << getDu(index,1,0,0,0.0,thesign,grid) << endl;
//   cout << tempval << " " << getDu(x,1,thesign,grid) << endl;
//   cout << tempval << " " << getu(index,0,0,0.0,thesign,grid) << endl;
//   cout << tempval << " " << getu(x,thesign,grid) << endl;

   double tempval2;
   tempval2 = c[N];
   for (r = 0; r < N; r++)
      tempval2 += c[r]*tempuval[r];
   cout << tempval2 << " " << getD2u(x,1,1,thesign,grid) << endl;
//   cout << tempval2 << " " << getDu(x,1,thesign,grid) << endl;
//   cout << tempval2 << " " << getu(x,thesign,grid) << endl;

   double theerr[N+1];
   tempval2 = c[N];
   tempval = c[N];
   theerr[N] = tempval-tempval2;
   cout << tempval << " " << tempval2 << " " << theerr[N] << endl;
   for (r = 0; r < N; r++)
   {
      tempval2 = c[r]*tempuval[r];
      if (evalarray(S,cindex[r]) < 0.0)
         tempval = c[r]*getu(cindex[r],0,0,0.0,-1,grid);
      else
         tempval = c[r]*getu(cindex[r],0,0,0.0,1,grid);
      theerr[r] = tempval-tempval2;
      cout << tempval << " " << tempval2 << " " << theerr[r] << endl;
   }

   free_matrix(B,N,N);
   free_matrix(LU,N,N);
   delete [] c;
   delete [] d;
}

void addtoheap(ZLHeapStruct &heap, int *index, double val)
{
   ZLHeapElt *parent, *child;
   int j, iparent, ichild, last;

   if (heap.tail != NULL)
      last = (*(heap.tail)).num+1;
   else
      last = 1;
   iparent = last>>1;
   ichild = last%2;
   if (iparent > 0)
      parent = heapgoto(heap,iparent);
   else
      parent = NULL;
   if (parent != NULL)
      if ((*parent).child[ichild] == NULL)
      {
         (*parent).child[ichild] = new ZLHeapElt;
         child = (*parent).child[ichild];
         (*child).parent = parent;
         (*child).child[0] = NULL;
         (*child).child[1] = NULL;
      }
      else
         child = (*parent).child[ichild];
   else
      if (heap.head == NULL)
      {
         heap.head = new ZLHeapElt;
         child = heap.head;
         (*child).parent = NULL;
         (*child).child[0] = NULL;
         (*child).child[1] = NULL;
      }
      else
         child = heap.head;

   (*child).num = last;
   for (j = 0; j < heap.dim; j++)
      (*child).index[j] = index[j];
   (*child).val = val;
   setvalarray(heap.mark,(*child).index,2);
   heap.tail = child;

   fixheapeltup(heap,child);
}

ZLHeapElt* heapgoto(ZLHeapStruct &heap, int num)
{
   ZLHeapElt *current;
   int n, k;

   if (num <= 0)
      return NULL;

   for (n = 1; n <= num; n*= 2);
   n /= 2;

   current = heap.head;
   k = num;
   while (n > 1 && current != NULL)
   {
      if (k >= n)
         k -= n;
      n /= 2;
      if (k >= n)
         current = (*current).child[1];
      else
         current = (*current).child[0];
   }

   return current;
}

void fixheapeltup(ZLHeapStruct &heap, ZLHeapElt *fix)
{
   ZLHeapElt *parent, *child, temp;
   char change = 1;
   int j;

   child = fix;
   parent = (*child).parent;
   while (change && parent != NULL)
   {
      if (fabs((*parent).val) > fabs((*child).val))
      {
         for (j = 0; j < heap.dim; j++)
         {
            temp.index[j] = (*parent).index[j];
            (*parent).index[j] = (*child).index[j];
            (*child).index[j] = temp.index[j];
            temp.val = (*parent).val;
            (*parent).val = (*child).val;
            (*child).val = temp.val;
         }
         child = parent;
         parent = (*parent).parent;
      }
      else
         change = 0;
   }
}

void fixheapeltdelete(ZLHeapStruct &heap, ZLHeapElt *del)
{
   ZLHeapElt *current, *current2;
   int j, temp;

   current = fixheapeltempty(heap,del);
   if (current != heap.tail)
   {
      for (j = 0; j < heap.dim; j++)
         (*current).index[j] = (*(heap.tail)).index[j];
      (*current).val = (*(heap.tail)).val;
      temp = (*(heap.tail)).num;
      current2 = (*(heap.tail)).parent;
      if (current2 != NULL)
      {
         for (j = 0; j < 2 && (*current2).child[j] != heap.tail; j++);
         (*current2).child[j] = NULL;
      }
      else
         heap.head = NULL;
      delete heap.tail;
      heap.tail = heapgoto(heap,temp-1);
      fixheapeltup(heap,current);
   }
   else
   {
      temp = (*(heap.tail)).num;
      current2 = (*(heap.tail)).parent;
      if (current2 != NULL)
      {
         for (j = 0; j < 2 && (*current2).child[j] != heap.tail; j++);
         (*current2).child[j] = NULL;
      }
      else
         heap.head = NULL;
      delete heap.tail;
      heap.tail = heapgoto(heap,temp-1);
   }
}

ZLHeapElt* fixheapeltempty(ZLHeapStruct &heap, ZLHeapElt *fix)
{
   ZLHeapElt *current, *child[2];
   int j, r;

   current = fix;
   child[0] = (*current).child[0];
   child[1] = (*current).child[1];
   while (child[0] != NULL && (*(child[0])).index[0] >= 0)
   {
      if (child[1] == NULL || (*(child[1])).index[0] < 0 ||
          fabs((*(child[0])).val) <= fabs((*(child[1])).val))
         r = 0;
      else
         r = 1;
      for (j = 0; j < heap.dim; j++)
         (*current).index[j] = (*(child[r])).index[j];
      (*current).val = (*(child[r])).val;
      current = child[r];
      child[0] = (*current).child[0];
      child[1] = (*current).child[1];
   }

   return current;
}

void readZLHeap(ZLHeapStruct heap)
{
   ZLHeapElt *current;
   int i;

   cout << "start read" << endl;
   for (i = 1,current = heapgoto(heap,i); current != NULL; 
        i++,current = heapgoto(heap,i))
       cout << "   " << (*current).num << " "
                     << (*current).index[0] << " "
                     << (*current).index[1] << " "
                     << (*current).index[2] << " "
                     << (*current).val << endl;
   cout << "end read" << endl;
}

void getnearestgrid(int **index, int &N, double *x, int maxpts, double maxdist, 
                    char ***tube, GridData &grid)
{
   ZLHeapStruct heap;
   int i, r, s, t;
   int rindex[grid.dim], sindex[grid.dim], tindex[grid.dim];
   double z[grid.dim];

   heap.head = NULL;
   heap.tail = NULL;
   heap.mark = tube;
   heap.dim = grid.dim;

   for (t = 0; t < grid.dim; t++)
      rindex[t] = (x[t]-grid.a[t])/grid.dx[t];

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= 1)
   {
      for (t = 0; t < grid.dim; t++)
      {
         sindex[t] = rindex[t]+tindex[t];
         if (sindex[t] < 0)
            sindex[t] = 0;
         else if (sindex[t] >= grid.nx[t])
            sindex[t] = grid.nx[t]-1;
      }
      if (evalarray(heap.mark,sindex) == 1)
      {
         sub2coord(z,sindex,grid);
         addtoheap(heap,sindex,getdist(z,x,grid.dim));
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > 1; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
//   readZLHeap(heap);
//   getchar();

   N = 0;
   while (heap.head != NULL && (*(heap.head)).index[0] >= 0 && 
          (maxpts < 0 || N < maxpts) && (maxdist < 0.0 || 
          (*(heap.head)).val <= maxdist))
   {
      for (r = 0; r < grid.dim; r++)
         index[N][r] = (*(heap.head)).index[r];
      fixheapeltdelete(heap,heap.head);
//      readZLHeap(heap);
//      getchar();
      for (t = 0; t < grid.dim; t++)
         sindex[t] = index[N][t];
      for (r = 0; r < grid.dim; r++)
         for (s = -1; s <= 1; s += 2)
         {
            sindex[r] = index[N][r]+s;
            if (sindex[r] >= 0 && sindex[r] <= grid.nx[r] && 
                evalarray(heap.mark,sindex) == 1)
            {
               sub2coord(z,sindex,grid);
               addtoheap(heap,sindex,getdist(z,x,grid.dim));
            }
            sindex[r] = index[N][r];
         }
//      readZLHeap(heap);
//      getchar();
      N++;
   }

   while (heap.tail != NULL)
   {
      setvalarray(heap.mark,(*(heap.tail)).index,1);
      fixheapeltdelete(heap,heap.tail);
   }
   for (t = 0; t < N; t++)
      setvalarray(heap.mark,index[t],1);
}

double getdist(double *z, double *x, int thedim)
{
   double val = 0.0;
   int r;

   for (r = 0; r < thedim; r++)
      val += (z[r]-x[r])*(z[r]-x[r]);
   val = sqrt(val);

   return val;
}

double bilinearinterp(double *x, double ***u, GridData &grid)
{
   int i, r, index[grid.dim], tindex[grid.dim];
   double value, temp, dindex[grid.dim];

   for (r = 0; r < grid.dim; r++)
      dindex[r] = (x[r]-grid.a[r])/grid.dx[r];
   for (r = 0; r < grid.dim; r++)
      index[r] = (int) dindex[r];

   value = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = index[i];
   while (tindex[0] <= index[0]+1)
   {
      temp = evalarray(u,tindex);
      for (r = 0; r < grid.dim; r++)
         temp *= 1.0-min(max(fabs(dindex[r]-tindex[r]),0.0),1.0);
      value += temp;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > index[i]+1; i--)
      {
         tindex[i] = index[i];
         (tindex[i-1])++;
      }
   }

   return value;
}

double weno6interpdirect(double *x, double ***u, GridData &grid)
{
   int r, s, i, j;
   int deg = 3;
   int index[grid.dim], tindex[grid.dim], rindex[grid.dim];
   double dindex[grid.dim];
   double tvalue, thesum, p, C, IS, alpha;
   double epsilon = 1.0e-6;
   double value[grid.dim][2*deg];

   for (r = 0; r < grid.dim; r++)
      dindex[r] = (x[r]-grid.a[r])/grid.dx[r];
   for (r = 0; r < grid.dim; r++)
      index[r] = (int) dindex[r];

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] < 2*deg)
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = min(max(index[r]-deg+1+tindex[r],0),grid.nx[r]);
      value[grid.dim-1][tindex[grid.dim-1]] = evalarray(u,rindex);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i >= 0 && tindex[i] >= 2*deg; i--)
      {
         j = deg-1;
         tvalue = 0.0;
         thesum = 0.0;
         for (s = 0; s < deg; s++)
         {
            p = value[i][s]+(dindex[i]-(index[i]-deg+1+s))*
                ((value[i][s+1]-value[i][s])+(dindex[i]-(index[i]-deg+2+s))*
                 ((value[i][s+2]-2.0*value[i][s+1]+value[i][s])/2.0+
                  (dindex[i]-(index[i]-deg+3+s))*
                  (value[i][s+3]-3.0*value[i][s+2]+3.0*value[i][s+1]-
                   value[i][s])/6.0));
            if (s == 0)
            {
               C = (index[i]+2-dindex[i])*(index[i]+3-dindex[i])/20.0;
               IS = (-3579.0*value[i][j+1]*value[i][j]+
                     2634.0*value[i][j+1]*value[i][j-1]-
                     683.0*value[i][j+1]*value[i][j-2]-
                     6927.0*value[i][j]*value[i][j-1]+
                     1854.0*value[i][j]*value[i][j-2]-
                     1659.0*value[i][j-1]*value[i][j-2]+
                     814.0*value[i][j+1]*value[i][j+1]+
                     4326.0*value[i][j]*value[i][j]+
                     2976.0*value[i][j-1]*value[i][j-1]+
                     244.0*value[i][j-2]*value[i][j-2])/180.0;
            }
            else if (s == 1)
            {
               C = (index[i]+3-dindex[i])*(dindex[i]-(index[i]-2))/10.0;
               IS = (-3777.0*value[i][j+1]*value[i][j]+
                     1074.0*value[i][j+1]*value[i][j-1]-
                     1269.0*value[i][j]*value[i][j-1]+
                     1986.0*value[i][j+1]*value[i][j+1]+
                     1986.0*value[i][j]*value[i][j]+
                     244.0*value[i][j-1]*value[i][j-1]+
                     244.0*value[i][j+2]*value[i][j+2]-
                     1269.0*value[i][j+2]*value[i][j+1]+
                     1074.0*value[i][j+2]*value[i][j]-
                     293.0*value[i][j+2]*value[i][j-1])/180.0;
            }
            else if (s == 2)
            {
               C = (dindex[i]-(index[i]-2))*(dindex[i]-(index[i]-1))/20.0;
               IS = (-3579.0*value[i][j+1]*value[i][j]+
                     4326.0*value[i][j+1]*value[i][j+1]+
                     814.0*value[i][j]*value[i][j]+
                     2976.0*value[i][j+2]*value[i][j+2]+
                     244.0*value[i][j+3]*value[i][j+3]-
                     683.0*value[i][j+3]*value[i][j]-
                     6927.0*value[i][j+2]*value[i][j+1]+
                     2634.0*value[i][j+2]*value[i][j]-
                     1659.0*value[i][j+3]*value[i][j+2]+
                     1854.0*value[i][j+3]*value[i][j+1])/180.0;
            }
            alpha = C/((epsilon+IS)*(epsilon+IS));
            thesum += alpha;
            tvalue += alpha*p;
         }
         tvalue /= thesum;
         if (i > 0)
         {
            value[i-1][tindex[i-1]] = tvalue;
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }

   return tvalue;
}

/* 
central difference approximate of gradient of S
grad[d][i][j][k], d = 1,2,3
not used
*/
void getallgrad(double ****grad, double ***S, GridData &grid)
{
   int i, r, s, tindex[grid.dim], sindex[grid.dim];
   double val;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      // sindex is nbr of tindex, if at boundary, sindex same as tindex
      for (r = 0; r < grid.dim; r++)
         sindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         s = 0;
         if (tindex[r]+1 <= grid.nx[r])
         {
            sindex[r] = tindex[r]+1;
            s++;
         }
         val = evalarray(S,sindex);
         if (tindex[r]-1 >= 0)
         {
            sindex[r] = tindex[r]-1;
            s++;
         }
         val -= evalarray(S,sindex);
         val /= s*grid.dx[r];
         setvalarray(grad[r],tindex,val);
         sindex[r] = tindex[r];
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}



int newtondir(double *x, double *x0, double *grad, double tol, double ***S, 
               double ****allgrad, GridData &grid)
{
   double f, fprime, thechange;
   int r, step;

   for (r = 0; r < grid.dim; r++)
      x[r] = x0[r];
   for (thechange = fabs(tol)+1.0,step = 0; fabs(thechange) > tol; step++)
   {
//      f = weno6interpdirect(x,S,grid);
      f = bilinearinterp(x,S,grid);
      fprime = 0.0;
      for (r = 0; r < grid.dim; r++)
//         fprime += weno6interpdirect(x,allgrad[r],grid)*grad[r];
         fprime += bilinearinterp(x,allgrad[r],grid)*grad[r];
      thechange = -f/fprime;
      for (r = 0; r < grid.dim; r++)
         x[r] += thechange*grad[r];
   }

   return step;
}

int regulafalsidir(double *x, double *x0, double *grad, double tol, double ***S, 
                   GridData &grid)
{
   double a, b, c, cold, fa, fb, fc, temp[grid.dim];
   int r, step, thesign;

   for (r = 0; r < grid.dim; r++)
      temp[r] = x0[r];

   a = 0.0;
   fa = weno6interpdirect(x0,S,grid);
   b  = 0.0;
   do 
   {
      b += 0.5*grid.mindx;
      for (r = 0; r < grid.dim; r++)
         x[r] = temp[r]+b*grad[r];
      fb = weno6interpdirect(x,S,grid);
   }
   while (fa*fb > 0.0);

   cold = b;
   for (c = a,step = 0; fabs(c-cold) > tol; step++)
   {
      cold = c;
      c = a-fa*(b-a)/(fb-fa);
      for (r = 0; r < grid.dim; r++)
         x[r] = temp[r]+c*grad[r];
      fc = weno6interpdirect(x,S,grid);
      if (globdebug)
         cout << "   approx at " << c << " with " << fc << endl;
      if ((fa > 0.0)+(fc > 0.0) == 1)
      {
         b = c;
         fb = fc;
      }
      else
      {
         a = c;
         fa = fc;
      }
   }

   return step;
}

/*
used in iim to find the nearest point on the interface
*/
void getnearestinterface(double *x, int *index, double ***S, GridData &grid) 
{
   int r, s, step, rindex[grid.dim];
   double grad[grid.dim];

   sub2coord(x,index,grid);
   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
   {
      grad[r] = 0.0;
      for (s = -1; s <= 1; s += 2)
      {
         rindex[r] = min(max(index[r]+s,0),grid.nx[r]);
         grad[r] += s*evalarray(S,rindex);
      }
      grad[r] /= 2.0*grid.dx[r];
      rindex[r] = index[r];
   }
   if (sqrt(getdotprod(grad,grad,grid.dim)) < 1.0e-11)
   {
      for (r = 0; r < grid.dim; r++)
         grad[r] = 0.0;
      grad[r] = 1.0;
   }
   if (evalarray(S,index) > 0.0)
      for (r = 0; r < grid.dim; r++)
         grad[r] = -grad[r];
      
//   step = newtondir(x,x,grad,1.0e-14,S,allgrad,grid);
   step = regulafalsidir(x,x,grad,1.0e-14,S,grid);
   if (globdebug)
   {
      cout << "took " << step << " number of root finding steps" << endl;
      cout << "value = " << bilinearinterp(x,S,grid) << " " 
                         << weno6interpdirect(x,S,grid) << endl;
      cout << "   at x = " << x[0] << " " << x[1] << " " << x[2] << endl;
      double y[grid.dim];
      sub2coord(y,index,grid);
      cout << "   which is dist = " << getdist(x,y,grid.dim) << " " << grid.dx[0] << endl;
      cout << "   with grad = " << grad[0] << " " << grad[1] << " " << grad[2] << endl;
      cout << "   and dist from origin = " << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << endl;
   }
}

double pythag(double a, double b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
   double absa, absb;
   absa = fabs(a);
   absb = fabs(b);
   if (absa > absb) 
      return absa*sqrt(1.0+(absb/absa)*(absb/absa));
   else 
      return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

void svdcmp(double **a, int m, int n, double w[], double **v)
{
   int flag, i, its, j, jj, k, l, nm;
   double anorm, c, f, g, h, s, scale, x, y, z, rv1[n];

/*
   cout << "   A = " << endl;
   for (i = 0; i < m; i++)
   {
      cout << "   ";
      for (j = 0; j < n; j++)
         cout << a[i][j] << " ";
      cout << endl;
   }
*/
   g = scale = anorm = 0.0; /* Householder reduction to bidiagonal form */
   for (i = 0; i < n; i++)
   {
      l = i+1;
      rv1[i] = scale*g;
      g = s = scale = 0.0;
      if (i < m)
      {
         for (k = i; k < m; k++)
            scale += fabs(a[k][i]);
         if (scale)
         {
            for (k = i; k < m; k++)
            {
               a[k][i] /= scale;
               s += a[k][i]*a[k][i];
            }
            f = a[i][i];
            g = -SIGN(sqrt(s),f);
            h = f*g-s;
            a[i][i] = f-g;
            for (j = l; j < n; j++)
            {
               for (s = 0.0,k = i; k < m; k++)
                  s += a[k][i]*a[k][j];
               f = s/h;
               for (k = i; k < m; k++)
                  a[k][j] += f*a[k][i];
            }
            for (k = i; k < m; k++)
               a[k][i] *= scale;
         }
      }
      w[i] = scale*g;
      g = s = scale = 0.0;
      if (i < m && i != n-1)
      {
         for (k = l; k < n; k++)
            scale += fabs(a[i][k]);
         if (scale)
         {
            for (k = l; k < n; k++)
            {
               a[i][k] /= scale;
               s += a[i][k]*a[i][k];
            }
            f = a[i][l];
            g = -SIGN(sqrt(s),f);
            h = f*g-s;
            a[i][l] = f-g;
            for (k = l; k < n; k++)
               rv1[k] = a[i][k]/h;
            for (j = l; j < m; j++)
            {
               for (s = 0.0,k = l; k < n; k++)
                  s += a[j][k]*a[i][k];
               for (k = l; k < n; k++)
                  a[j][k] += s*rv1[k];
            }
            for (k = l; k < n; k++)
               a[i][k] *= scale;
         }
      }
      anorm = DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
   }
   for (i = n-1; i >= 0; i--) /* Accumulation of right-hand transformations. */
   {
      if (i < n-1)
      {
         if (g)
         {
            for (j = l; j < n; j++) /* Double division to avoid possible underflow. */
               v[j][i] = (a[i][j]/a[i][l])/g;
            for (j = l; j < n; j++)
            {
               for (s = 0.0,k = l; k < n; k++)
                  s += a[i][k]*v[k][j];
               for (k = l; k < n; k++)
                  v[k][j] += s*v[k][i];
            }
         }
         for (j = l; j < n; j++)
            v[i][j] = v[j][i] = 0.0;
      }
      v[i][i] = 1.0;
      g = rv1[i];
      l = i;
   }
   for (i = IMIN(m-1,n-1); i >= 0; i--) /* Accumulation of left-hand transformations. */
   {
      l = i+1;
      g = w[i];
      for (j = l; j < n; j++)
         a[i][j] = 0.0;
      if (g)
      {
         g = 1.0/g;
         for (j = l; j < n; j++)
         {
            for (s = 0.0,k = l; k < m; k++)
               s += a[k][i]*a[k][j];
            f = (s/a[i][i])*g;
            for (k = i; k < m; k++)
               a[k][j] += f*a[k][i];
         }
         for (j = i; j < m; j++)
            a[j][i] *= g;
      }
      else for (j = i; j < m; j++)
         a[j][i]=0.0;
      ++a[i][i];
   }
   for (k = n-1; k >= 0; k--) /* Diagonalization of the bidiagonal form. */
   {
      for (its = 1; its <= 30; its++)
      {
         flag = 1;
         for (l = k; l >= 0; l--) /* Test for splitting. */
         {
            nm = l-1; /* Note that rv1[1] is always zero. */
            if ((double)(fabs(rv1[l])+anorm) == anorm)
            {
               flag = 0;
               break;
            }
            if ((double)(fabs(w[nm])+anorm) == anorm)
               break;
         }
         if (flag)
         {
            c=0.0; /* Cancellation of rv1[l], if l > 1. */
            s=1.0;
            for (i = l; i <= k; i++)
            {
               f = s*rv1[i];
               rv1[i] = c*rv1[i];
               if ((double)(fabs(f)+anorm) == anorm)
                  break;
               g = w[i];
               h = pythag(f,g);
               w[i] = h;
               h = 1.0/h;
               c = g*h;
               s = -f*h;
               for (j = 0; j < m; j++)
               {
                  y = a[j][nm];
                  z = a[j][i];
                  a[j][nm] = y*c+z*s;
                  a[j][i] = z*c-y*s;
               }
            }
         }
         z = w[k];
         if (l == k) /* Convergence. */
         {
            if (z < 0.0) /* Singular value is made nonnegative. */
            {
               w[k] = -z;
               for (j = 0; j < n; j++)
                  v[j][k] = -v[j][k];
            }
            break;
         }
         if (its == 30)
            printf("no convergence in 30 svdcmp iterations");
         x = w[l]; /* Shift from bottom 2-by-2 minor. */
         nm = k-1;
         y = w[nm];
         g = rv1[nm];
         h = rv1[k];
         f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
         g = pythag(f,1.0);
         f = ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
         c = s = 1.0; /* Next QR transformation: */
         for (j = l; j <= nm; j++)
         {
            i = j+1;
            g = rv1[i];
            y = w[i];
            h = s*g;
            g = c*g;
            z = pythag(f,h);
            rv1[j] = z;
            c = f/z;
            s = h/z;
            f = x*c+g*s;
            g = g*c-x*s;
            h = y*s;
            y *= c;
            for (jj = 0; jj < n; jj++)
            {
               x = v[jj][j];
               z = v[jj][i];
               v[jj][j] = x*c+z*s;
               v[jj][i] = z*c-x*s;
            }
            z = pythag(f,h);
            w[j] = z; /* Rotation can be arbitrary if z = 0. */
            if (z)
            {
               z = 1.0/z;
               c = f*z;
               s = h*z;
            }
            f = c*g+s*y;
            x = c*y-s*g;
            for (jj = 0; jj < m; jj++)
            {
               y = a[jj][j];
               z = a[jj][i];
               a[jj][j] = y*c+z*s;
               a[jj][i] = z*c-y*s;
            }
         }
         rv1[l] = 0.0;
         rv1[k] = f;
         w[k] = x;
      }
   }

/*
   cout << "   U = " << endl;
   for (i = 0; i < m; i++)
   {
      cout << "   ";
      for (j = 0; j < n; j++)
         cout << a[i][j] << " ";
      cout << endl;
   }
   cout << "   S = " << endl;
   cout << "   ";
   for (i = 0; i < n; i++)
      cout << w[i] << " ";
   cout << endl;
   cout << "   V = " << endl;
   for (i = 0; i < n; i++)
   {
      cout << "   ";
      for (j = 0; j < n; j++)
         cout << v[i][j] << " ";
      cout << endl;
   }
   double **B = matrix(m-1,n-1);
   for (i = 0; i < m; i++)
   {
      for (j = 0; j < n; j++)
      {
         B[i][j] = 0.0;
         for (k = 0; k < n; k++)
            B[i][j] += a[i][k]*w[k]*v[j][k];
      }
   }
   cout << "   A = " << endl;
   for (i = 0; i < m; i++)
   {
      cout << "   ";
      for (j = 0; j < n; j++)
         cout << B[i][j] << " ";
      cout << endl;
   }
   free_matrix(B,m-1,n-1);
*/
}

double svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[],
              double thresh)
{
   int jj, j, i;
   double s, tsh, eps, tmp[n];

   eps = numeric_limits<double>::epsilon();
   tsh = (thresh >= 0.0 ? thresh : 0.5*sqrt(m+n+1.0)*w[0]*eps);
//   tsh = (thresh >= 0.0 ? thresh : 10.0*eps);
   for (j = 0; j < n; j++)
   {
      s = 0.0;
      if (w[j] > tsh)
      {
         for (i = 0; i < m; i++)
            s += u[i][j]*b[i];
         s /= w[j];
      }
      tmp[j] = s;
   }
   for (j = 0; j < n; j++)
   {
      s = 0.0;
      for (jj = 0; jj < n; jj++)
         s += v[j][jj]*tmp[jj];
      x[j] = s;
   }

   double residual[m], tresidual;
   int p, q;

   tresidual = 0.0;
   for (i = 0; i < m; i++)
   {
      residual[i] = b[i];
      for (p = 0; p < n; p++)
         for (q = 0; q < n; q++)
            residual[i] -= u[i][p]*w[p]*v[q][p]*x[q];
      tresidual += residual[i]*residual[i];
   }
   tresidual = sqrt(tresidual);
   if (globdebug)
      cout << "residual length = " << tresidual << endl;
   return tresidual;
}

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

void testZL(char yesC, int add, double ***S, char ***tube, PBData &pb, 
            GridData &grid)
{
   int M = 2+grid.dim+(grid.dim*(grid.dim+1))/2;
   int N = 1;
   int i, j, k, r, s, t, rstar, sstar;
   for (r = 0; r < grid.dim; r++)
      N *= 3;
   int L;
   double **B, *d, *c, **LU, *w;
   int *PLR, *PLC;
   int index[grid.dim], rindex[grid.dim], **cindex; 
   double up0, upum, uxp0[grid.dim], **uxpuxm = matrix(grid.dim-1,grid.dim-1),
          **uxxp0 = matrix(grid.dim-1,grid.dim-1), 
          ***uxxpuxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1),
          ****uxxpuxxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1,grid.dim-1);
   double ehere, ethere, thesign, alpha, tempalpha; 
   double x[grid.dim], normal[grid.dim], tempnormal[grid.dim], junk[grid.dim];
   int two2one[grid.dim][grid.dim];
   double theerr[M];
   for (r = 0; r < M; r++)
      theerr[r] = 0.0;
   int theerri[M][grid.dim];
   for (r = 0; r < M; r++)
      for (s = 0; s < grid.dim; s++)
         theerri[r][s] = -1;
   double theerrx[M][grid.dim];
   for (r = 0; r < M; r++)
      for (s = 0; s < grid.dim; s++)
         theerrx[r][s] = -2.0;

   for (i = 0; i < grid.dim; i++)
      index[i] = 0;
   while (index[0] <= grid.nx[0])
   {
      int thestatus = 1;
      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];
      for (r = 0; r < grid.dim && thestatus == 1; r++)
      {
         for (s = -1; s <= 1 && thestatus == 1; s += 2)
         {
            rindex[r] = index[r]+s;
            if (rindex[r] >= 0 && rindex[r] <= grid.nx[r])
               if (evalarray(S,index)*evalarray(S,rindex) < 0.0)
                  thestatus = 2;
               else;
            else
               thestatus = 0;
         }
         rindex[r] = index[r];
      }

      if (thestatus == 2)
      {
//         cout << index[0] << " " << index[1] << " " << index[2] << endl;
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
            two2one[r][r] = r;
         two2one[0][1] = grid.dim;
         two2one[1][2] = grid.dim+1;
         two2one[0][2] = grid.dim+2;
         for (s = 0; s < grid.dim; s++)
            for (r = 0; r < s; r++)
               two2one[s][r] = two2one[r][s];
     
         getnearestinterface(x,index,S,grid);
         getiimjumps(up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm,x,thesign,S,pb,grid);
         if (yesC)
         {
            cindex = imatrix(M-1+add-1,grid.dim-1);
            getnearestgrid(cindex,N,x,M-1+add,-1.0,tube,grid);
            L = N+1;
         }
         else
         {
            cindex = imatrix(M+add-1,grid.dim-1);
            getnearestgrid(cindex,N,x,M+add,-1.0,tube,grid);
            L = N;
         }

//         if (index[0] == 29 && index[1] == 20 && index[2] == 20)
//            for (r = 0; r < L; r++)
//            {
//               for (s = 0; s < grid.dim; s++)
//                  cout << cindex[r][s] << " ";
//               cout << endl;
//            }

         B = matrix(L-1,L-1); 
         d = new double[L]; 
         c = new double[L];
         LU = matrix(L-1,L-1);
         PLR = new int[L]; 
         PLC = new int[L];
         w = new double[L];

         getiimstencilmx(B,M,N,x,thesign,cindex,yesC,up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,
                         uxxpuxxm,S,pb,grid);
//         getiimgridstencilmx(B,M,N,x,index,thesign,cindex,yesC,up0,upum,uxp0,uxpuxm,uxxp0,
//                             uxxpuxm,uxxpuxxm,S,pb,grid);
         for (r = 0; r < M; r++)
            d[r] = 0.0;
         for (k = 0; k < M; k++)
         {
            d[k] = 1.0;
            
            double **AA = matrix(M-1,L-1), dd[M];
            for (r = 0; r < M; r++)
            {
               for (s = 0; s < L; s++)
                  AA[r][s] = B[r][s];
               dd[r] = d[r];
            }
//            double smallnum = 1.0;
//            t = 0;
//            for (s = 0; s < N && t != grid.dim; s++)
//               for (t = 0; t < grid.dim && cindex[s][t] == index[t]; t++);
//            s--;
//            for (r = 0; r < M; r++)
//               AA[r][s] /= smallnum;
      
            svdcmp(AA,M,L,w,LU);
            svbksb(AA,w,LU,M,L,dd,c,1.0e-14);
//            c[s] /= smallnum;
            free_matrix(AA,M-1,L-1);
     
//            if (index[0] == 29 && index[1] == 20 && index[2] == 20)
//            {
//               for (r = 0; r < N; r++)
//                  cout << c[r] << " ";
//               cout << endl;
//            }
   
            double tempval, exactval;
   
            if (yesC)
               tempval = c[N];
            else
               tempval = 0.0;
            for (r = 0; r < N; r++)
               if (evalarray(S,cindex[r]) < 0.0)
                  tempval += c[r]*getu(cindex[r],0,0,0.0,-1,grid);
               else
                  tempval += c[r]*getu(cindex[r],0,0,0.0,1,grid);
            if (k == 0)
            {
               exactval = 1.0;
//               cout << "error for const" << endl;
            }
            else if (k == 1)
            {
               exactval = getu(x,thesign,grid);
//               cout << "error for u" << endl;
            }
            else if (k < 2+grid.dim)
            {
               exactval = getDu(x,k-2,thesign,grid);
//               cout << "error for Du " << k-2 << endl;
            }
            else
            {
               for (r = 0; r < grid.dim; r++)
                  for (s = r; s < grid.dim; s++)
                     if (two2one[r][s] == k-2-grid.dim)
                     {
                        exactval = getD2u(x,r,s,thesign,grid);
//                        cout << "error for D2u " << r << " " << s << endl;
                     }
            }
            if (index[0] == 29 && index[1] == 20 && index[2] == 20)
               cout << "   " << tempval << " " << exactval << " " 
                    << fabs(tempval-exactval) << endl;
            if (fabs(tempval-exactval) > theerr[k])
            {
               theerr[k] = fabs(tempval-exactval);
               for (r = 0; r < grid.dim; r++)
                  theerri[k][r] = index[r];
               for (r = 0; r < grid.dim; r++)
                  theerrx[k][r] = x[r];
               if (fabs(tempval-exactval) > 4.5)
                  cout << x[0] << " " << x[1] << " " << x[2] << " " 
                       << fabs(tempval-exactval) << endl;
            }
            d[k] = 0.0;
         }

         free_matrix(B,L-1,L-1);
         free_matrix(LU,L-1,L-1);
         free_matrix(cindex,N-1,grid.dim-1);
         delete [] c;
         delete [] d;
         delete [] PLR;
         delete [] PLC;
         delete [] w;

//         getchar();
      }

      (index[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && index[i] > grid.nx[i]; i--)
      {
         index[i] = 0;
         (index[i-1])++;
      }
   }

   cout << "Error from testZL is" << endl << "   ";
   for (r = 0; r < M; r++)
      cout << theerr[r] << " ";
   cout << endl;
   for (s = 0; s < grid.dim; s++)
   {
      for (r = 0; r < M; r++)
         cout << theerri[r][s] << " ";
      cout << endl;
   }
   for (r = 0; r < M; r++)
   {
      for (s = 0; s < grid.dim; s++)
         index[s] = theerri[r][s];
      if (evalarray(S,index) < 0.0)
         cout << "-1 ";
      else
         cout << "1 ";
   }
   cout << endl;
   for (s = 0; s < grid.dim; s++)
   {
      for (r = 0; r < M; r++)
         printf("%4.16f ", theerrx[r][s]);
//         cout << theerrx[r][s] << " ";
      cout << endl;
   }

/*
   x[0] = 0.4999999999999899;
   x[1] = 0.0;
   x[2] = 0.0;
   thesign = -1.0;
//   testZLatx(theerr,x,thesign,yesC,add,S,tube,pb,grid);
   testZLatx(theerr,x,thesign,yesC,10,S,tube,pb,grid);
   cout << "theerr again = ";
   for (r = 0; r < M; r++)
      cout << theerr[r] << " ";
   cout << endl;
   exit(1);
*/
}

void testZLatx(double *theerr, double *x, double thesign, char yesC, int add, 
               double ***S, char ***tube, PBData &pb, GridData &grid)
{
   int M = 2+grid.dim+(grid.dim*(grid.dim+1))/2;
   int N = 1;
   int i, j, k, r, s, t, rstar, sstar;
   for (r = 0; r < grid.dim; r++)
      N *= 3;
   int L;
   double **B, *d, *c, **LU, *w;
   int *PLR, *PLC;
   int **cindex; 
   double up0, upum, uxp0[grid.dim], **uxpuxm = matrix(grid.dim-1,grid.dim-1),
          **uxxp0 = matrix(grid.dim-1,grid.dim-1), 
          ***uxxpuxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1),
          ****uxxpuxxm = matrix(grid.dim-1,grid.dim-1,grid.dim-1,grid.dim-1);
   double ehere, ethere, alpha, tempalpha; 
   double normal[grid.dim], tempnormal[grid.dim], junk[grid.dim];
   int two2one[grid.dim][grid.dim];

   for (r = 0; r < grid.dim; r++)
      two2one[r][r] = r;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (s = 0; s < grid.dim; s++)
      for (r = 0; r < s; r++)
         two2one[s][r] = two2one[r][s];
     
   getiimjumps(up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,uxxpuxxm,x,thesign,S,pb,grid);
   if (yesC)
   {
      cindex = imatrix(M-1+add-1,grid.dim-1);
      getnearestgrid(cindex,N,x,M-1+add,-1.0,tube,grid);
      L = N+1;
   }
   else
   {
      cindex = imatrix(M+add-1,grid.dim-1);
      getnearestgrid(cindex,N,x,M+add,-1.0,tube,grid);
      L = N;
   }

//   for (r = 0; r < L; r++)
//   {
//      for (s = 0; s < grid.dim; s++)
//         cout << cindex[r][s] << " ";
//      cout << endl;
//   }

   B = matrix(L-1,L-1); 
   d = new double[L]; 
   c = new double[L];
   LU = matrix(L-1,L-1);
   PLR = new int[L]; 
   PLC = new int[L];
   w = new double[L];

   getiimstencilmx(B,M,N,x,thesign,cindex,yesC,up0,upum,uxp0,uxpuxm,uxxp0,uxxpuxm,
                   uxxpuxxm,S,pb,grid);
   for (r = 0; r < M; r++)
      d[r] = 0.0;
   for (k = 0; k < M; k++)
   {
      d[k] = 1.0;
      
      double **AA = matrix(M-1,L-1), dd[M];
      for (r = 0; r < M; r++)
      {
         for (s = 0; s < L; s++)
            AA[r][s] = B[r][s];
         dd[r] = d[r];
      }
      double smallnum = 1.0e7;
      for (s = 0; s < N; s++)
         if (evalarray(S,cindex[s]) < 0.0) 
            for (r = 0; r < M; r++)
               AA[r][s] /= smallnum;

      svdcmp(AA,M,L,w,LU);
      svbksb(AA,w,LU,M,L,dd,c,1.0e-14);
      for (s = 0; s < N; s++)
         if (evalarray(S,cindex[s]) < 0.0) 
            c[s] /= smallnum;
      free_matrix(AA,M-1,L-1);

//      for (r = 0; r < N; r++)
//         cout << c[r] << " ";
//      cout << endl;

      double tempval, exactval;

      if (yesC)
         tempval = c[N];
      else
         tempval = 0.0;
      for (r = 0; r < N; r++)
         if (evalarray(S,cindex[r]) < 0.0)
            tempval += c[r]*getu(cindex[r],0,0,0.0,-1,grid);
         else
            tempval += c[r]*getu(cindex[r],0,0,0.0,1,grid);
      if (k == 0)
      {
         exactval = 1.0;
//         cout << "error for const" << endl;
      }
      else if (k == 1)
      {
         exactval = getu(x,thesign,grid);
//         cout << "error for u" << endl;
      }
      else if (k < 2+grid.dim)
      {
         exactval = getDu(x,k-2,thesign,grid);
//         cout << "error for Du " << k-2 << endl;
      }
      else
      {
         for (r = 0; r < grid.dim; r++)
            for (s = r; s < grid.dim; s++)
               if (two2one[r][s] == k-2-grid.dim)
               {
                  exactval = getD2u(x,r,s,thesign,grid);
//                  cout << "error for D2u " << r << " " << s << endl;
               }
      }
      cout << "   " << tempval << " " << exactval << " " 
           << fabs(tempval-exactval) << endl;
      theerr[k] = fabs(tempval-exactval);
      d[k] = 0.0;
      if (k == 5)
      {
         cout << "c for uxx:" << endl;
         for (r = 0; r < N; r++)
            cout << c[r] << " ";
         cout << endl;

         double terr[N], y[grid.dim], value;
         cout << "taylor error for uxx:" << endl;
         for (r = 0; r < N; r++)
         {
            sub2coord(y,cindex[r],grid);
            if (evalarray(S,cindex[r])*thesign >= 0.0) 
            {
               terr[r] = getu(x,thesign,grid);
               for (s = 0; s < grid.dim; s++)
                  terr[r] += (y[s]-x[s])*getDu(x,s,thesign,grid);
               for (s = 0; s < grid.dim; s++)
                  for (t = 0; t < grid.dim; t++)
                     terr[r] += 0.5*(y[s]-x[s])*(y[t]-x[t])*getD2u(x,s,t,thesign,grid);
               cout << "same " << terr[r] << " " << getu(y,thesign,grid) << " " 
                    << terr[r]-getu(y,thesign,grid) << " " << thesign << endl;
               terr[r] = terr[r]-getu(y,thesign,grid);
            }
            else
            {
/*
               terr[r] = getu(x,-thesign,grid);
               for (s = 0; s < grid.dim; s++)
                  terr[r] += (y[s]-x[s])*getDu(x,s,-thesign,grid);
               for (s = 0; s < grid.dim; s++)
                  for (t = 0; t < grid.dim; t++)
                     terr[r] += 0.5*(y[s]-x[s])*(y[t]-x[t])*getD2u(x,s,t,-thesign,grid);
               cout << "other side first " << terr[r] << " " << getu(y,-thesign,grid) 
                    << " " << terr[r]-getu(y,-thesign,grid) << endl;
*/

               value = upum*getu(x,thesign,grid)+up0;
               terr[r] = value;
               for (s = 0; s < grid.dim; s++)
               {
                  value = uxp0[s];
                  for (t = 0; t < grid.dim; t++)
                     value += uxpuxm[s][t]*getDu(x,t,thesign,grid);
                  terr[r] += (y[s]-x[s])*value;
               }
               for (s = 0; s < grid.dim; s++)
                  for (t = 0; t < grid.dim; t++)
                  {
                     value = uxxp0[s][t];
                     for (j = 0; j < grid.dim; j++)
                        value += uxxpuxm[s][t][j]*getDu(x,j,thesign,grid);
                     for (int m = 0; m < grid.dim; m++)
                        for (j = 0; j < grid.dim; j++)
                           value += uxxpuxxm[s][t][m][j]*getD2u(x,m,j,thesign,grid);
                     terr[r] += 0.5*(y[s]-x[s])*(y[t]-x[t])*value;
                  }
               cout << "other " << terr[r] << " " << getu(y,-thesign,grid) << " " 
                    << terr[r]-getu(y,-thesign,grid) << " " << -thesign << endl;
               terr[r] = terr[r]-getu(y,-thesign,grid);
            }
         }
         value = 0.0;
         for (r = 0; r < N; r++)
            value += terr[r]*c[r];
         cout << "error should be " << value << endl;
      }
   }

   free_matrix(B,L-1,L-1);
   free_matrix(LU,L-1,L-1);
   free_matrix(cindex,N-1,grid.dim-1);
   delete [] c;
   delete [] d;
   delete [] PLR;
   delete [] PLC;
   delete [] w;

//   getchar();
}
