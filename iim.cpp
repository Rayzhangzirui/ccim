#include "iim.h"
#include "numerics.h"
#include "interface.h"

#include <iostream>
#include <cmath>
extern int globdebug;
using namespace std;

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
            rindex[r] = fmin(fmax(sindex[r]+s,0),grid.nx[r]);
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
                  rindex[r] = fmin(fmax(sindex[r]+s,0),grid.nx[r]);
                  val += evalarray(S,rindex)-evalarray(S,sindex);
               }
               val /= grid.dx[r]*grid.dx[r];
            }
            else
            {
               for (s = -1; s <= 1; s += 2)
               {
                  rindex[r] = fmin(fmax(sindex[r]+s,0),grid.nx[r]);
                  for (t = -1; t <= 1; t += 2)
                  {
                     rindex[m] = fmin(fmax(sindex[m]+t,0),grid.nx[r]);
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
