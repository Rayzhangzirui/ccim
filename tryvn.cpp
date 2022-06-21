#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <string>
#include <limits>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cassert>


#include "tryvn.h"
#include "advance.h"
#include "helper.h"
#include "icim.h"
#include "interface.h"
#include "iim.h"
#include "gmres.h"
#include "cim12.h"
#include "cim345cond.h"
#include "global.h"



using namespace std;


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
          // if index is boundary. 
         // This never happens. interiorptsmall is called only at interior point in linearsys
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
        checkanswer(pb.psi,u,grid);
      }
      if (globcim == 5)
      {
         checkcim5Du(pb.psi,u,pb,grid);
         checkallDu(pb.psi,u,pb,grid);
      }
      else if (globcim == 345 || globcim == 0 ||globcim == 6)
      {
         if (globsmall)
            checkDuStorage(pb.psi,Dusmall,smallsize,u,pb,grid);
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
          checkDuStorage(pb.psi,Dusmall,smallsize,u,pb,grid);
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
      checkanswer(pb.psi,u,grid);
      if (globcim == 5)
      {
         checkcim5Du(pb.psi,u,pb,grid);
         checkallDu(pb.psi,u,pb,grid);
      }
      else if (globcim == 345 || globcim == 0 || globcim == 6)
      {
         if (globsmall)
            checkDuStorage(pb.psi,Dusmall,smallsize,u,pb,grid);
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

