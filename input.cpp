#include "input.h"
#include "global.h"
#include "helper.h"
#include <cmath>
#include <iostream>
#include <sstream>

using namespace std;

#include "extratest.h"
#include "customtest.h"

double DBC(double *x, int thedim, double thetime)
{
// usually assume boundary uses outside(plus side) 
   if (globtestnum == 0)
   {
      double value = 0.0;
      
      for (int j = 0; j < thedim; j++)
         value += x[j]*x[j];
      value = -1.0/(1.0+value);

      return value;
   }
   else if (globtestnum == 1)
   {
      
      double value = 0.0;
      
      for (int j = 0; j < thedim; j++)
         value += x[j]*x[j];
   
      return value;
   }
   else if (globtestnum == 2)
   {
      
      double value = 0.0;
      
      for (int j = 0; j < thedim; j++)
         value += x[j]*x[j];
   
      return exp(value);
   }
   else if (globtestnum == 3 || globtestnum == 4)
   {
      
      double value = 0.0, epsp = EPSILONP, epsm = EPSILONM;
      
      for (int j = 0; j < thedim; j++)
         value += x[j]*x[j];
   
      return sqrt(value+fabs(1.0-epsp/epsm));
   }
   else if(globtestnum == 12 ||globtestnum == 11 || globtestnum == 10){
      return getuTest(x,1);
   }
   else if(globtestnum == 9){
      return dbc_custom(x,thedim,0);
   }

   return 0.0;
}



double getf(double *x, double thesign, PBData &pb, GridData &grid){
  double value = 0.0;
  double radius2 = 0.0;
  for (int r = 0; r < grid.dim; r++)
    radius2 += x[r]*x[r];

   if (globtestnum == 0)
   {
      
      if (thesign < 0.0)
      {
         for (int r = 0; r < grid.dim; r++)
            value -= pb.epsilonm*
                     (8.0*x[r]*x[r]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2))-
                      2.0/((1.0+radius2)*(1.0+radius2)));
      }
      else
      {
         for (int r = 0; r < grid.dim; r++)
            value -= pb.epsilonp*
                     (-8.0*x[r]*x[r]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2))+
                      2.0/((1.0+radius2)*(1.0+radius2)));
      }
   }
   else if (globtestnum == 1)
   {
      
      value = -2.0*grid.dim*pb.epsilonp;
   }
   else if (globtestnum == 2)
   {
      
      if (thesign < 0.0)
         value = -pb.epsilonm*(2.0*grid.dim*pb.epsilonp/pb.epsilonm+
                               4.0*(pb.epsilonp/pb.epsilonm)*(pb.epsilonp/pb.epsilonm)*
                               radius2)*exp(pb.epsilonp/pb.epsilonm*radius2-
                                            (pb.epsilonp/pb.epsilonm-1.0)*
                                            grid.radius*grid.radius);
      else
         value = -pb.epsilonp*(2.0*grid.dim+4.0*radius2)*exp(radius2);

   }
   else if (globtestnum == 3)
   {
      double A = 1.0, B = fabs(1.0-pb.epsilonp/pb.epsilonm), 
             C = pb.epsilonp/pb.epsilonm, 
             D = grid.radius*grid.radius*(1.0-pb.epsilonp/pb.epsilonm)+B;

      if (thesign < 0.0)
         value = -pb.epsilonm*C*(grid.dim-C*radius2/(C*radius2+D))/sqrt(C*radius2+D);
      else
         value = -pb.epsilonp*A*(grid.dim-A*radius2/(A*radius2+B))/sqrt(A*radius2+B);

   }
   else if (globtestnum == 4)
   {
      
      double A = 1.0, B = fabs(1.0-pb.epsilonp/pb.epsilonm), 
             C = pb.epsilonp/pb.epsilonm, 
             D = grid.radius*grid.radius*(1.0-pb.epsilonp/pb.epsilonm)+B;

      if (thesign < 0.0)
         value = -pb.epsilonm*C*(grid.dim-C*radius2/(C*radius2+D))/sqrt(C*radius2+D)+
                 pb.epsilonm*sin(x[0])*sqrt(C*radius2+D);
      else
         value = -pb.epsilonp*A*(grid.dim-A*radius2/(A*radius2+B))/sqrt(A*radius2+B)+
                 pb.epsilonp*cos(x[grid.dim-1])*sqrt(A*radius2+B);

      
   }
   else if (globtestnum == 12 ||globtestnum == 11 || globtestnum == 10)
   {  
      value = getfTest(x,thesign,pb);
   }
   else if (globtestnum == 9)
   {  
      value = getf_custom(x,thesign,pb, grid);
   }
   else{
      cerr<<"undefined globtestnum in getf"<<endl;
      exit(1);
   }


   // with a term
  if(globtesta == 0){
    // without a term
     value+=0.0;
  }else{
    value += geta(x, thesign, pb, grid) * getu(x,thesign, grid);
  }
    
    
  return value;
}


// getf between gridpoint use alpha
double getf(int *index, int rstar, int sstar, double alpha, double thesign, PBData &pb, 
            GridData &grid)
{
    double x[grid.dim];
    sub2coord(x,index,grid);
    x[rstar] += sstar*alpha*grid.dx[rstar];
    return getf(x, thesign, pb, grid);
}

double getu(double* x, double thesign, GridData &grid)
{
   double value = 0.0;     
   if (globtestnum == 0)
   {
      double radius2 = 0.0;
      for (int r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
      if (thesign < 0.0)
         value = 1.0/(1.0+radius2);
      else
         value = -1.0/(1.0+radius2);

      return value;
   }
   else if (globtestnum == 1)
   {
      double epsp = EPSILONP, epsm = EPSILONM;
      double radius2 = 0.0;
      for (int r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
      if (thesign < 0.0)
         value = epsp/epsm*radius2-(epsp/epsm-1.0)*pow(grid.radius0,2)*exp(4.0*(1.0-epsp/epsm)*grid.t);
      else
         value = radius2;

      return value;
   }
   else if (globtestnum == 2)
   {

      double epsp = EPSILONP, epsm = EPSILONM;

      double radius2 = 0.0;
      for (int r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
      if (thesign < 0.0)
         value = exp(epsp/epsm*radius2-(epsp/epsm-1.0)*grid.radius*grid.radius);
      else
         value = exp(radius2);

      return value;
   }
   else if (globtestnum == 3 || globtestnum == 4)
   {
      
      double epsp = EPSILONP, epsm = EPSILONM;
      double A = 1.0, B = fabs(1.0-epsp/epsm), C = epsp/epsm, 
             D = grid.radius*grid.radius*(1.0-epsp/epsm)+B;
      double radius2 = 0.0;
      for (int r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
      if (thesign < 0.0)
         value = sqrt(C*radius2+D);
      else
         value = sqrt(A*radius2+B);

      return value;
   }
   else if (globtestnum == 12 ||globtestnum == 11 || globtestnum == 10)
   {  
      return getuTest(x,thesign);
   
   }else{
      cerr<<"undefined globtestnum in getu"<<endl;
      exit(1);
   }
   return value;
}



double getu(int *index, int rstar, int sstar, double alpha, double thesign, GridData &grid)
{
    double x[grid.dim];
    sub2coord(x,index,grid);
    x[rstar] += sstar*alpha*grid.dx[rstar];
    return getu(x, thesign, grid);
}


double getDu(double *x, int s, double thesign, GridData &grid)
{
   if (globtestnum == 0)
   {
      int r;
      double value, radius2;

      radius2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
      if (thesign < 0.0)
      {
         value = -2.0*x[s]/((1.0+radius2)*(1.0+radius2));
      }
      else
      {
         value = 2.0*x[s]/((1.0+radius2)*(1.0+radius2));
      }
   
      return value;
   }
   else if (globtestnum == 1)
   {
      int r;
      double value;
      double epsp = EPSILONP, epsm = EPSILONM;

      if (thesign < 0.0)
         value = 2.0*epsp/epsm*x[s];
      else
         value = 2.0*x[s];
   
      return value;
   }
   else if (globtestnum == 2)
   {
      int r;
      double value, radius2;
      double epsp = EPSILONP, epsm = EPSILONM;

      radius2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];

      if (thesign < 0.0)
         value = 2.0*epsp/epsm*x[s]*exp(epsp/epsm*radius2-
                                        (epsp/epsm-1.0)*grid.radius*grid.radius);
      else
         value = 2.0*x[s]*exp(radius2);
   
      return value;
   }
   else if (globtestnum == 3 || globtestnum == 4)
   {
      int r;
      double value, radius2;
      double epsp = EPSILONP, epsm = EPSILONM;
      double A = 1.0, B = fabs(1.0-epsp/epsm), 
             C = epsp/epsm, D = grid.radius*grid.radius*(1.0-epsp/epsm)+B;

      radius2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];

      if (thesign < 0.0)
         value = C*x[s]/sqrt(C*radius2+D);
      else
         value = A*x[s]/sqrt(A*radius2+B);
   
      return value;
   }
   else if (globtestnum == 12 ||globtestnum == 11 || globtestnum == 10)
   {  
      return getDuTest(x,s,thesign);
   }


   cerr << "undefined Du" << endl;
   exit(-1);
   return 0.0;
}

double getDu(int *index, int s, int rstar, int sstar, double alpha, double thesign, 
             GridData &grid)
{
   double x[grid.dim];
   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];
   return getDu(x, s, thesign,grid);
}

double getD2u(double *x, int r, int s, double thesign, GridData &grid)
{
   if (globtestnum == 0)
   {
      int t;
      double value, radius2;

      radius2 = 0.0;
      for (t = 0; t < grid.dim; t++)
         radius2 += x[t]*x[t];
      if (thesign < 0.0)
         if (r == s)
            value = 8.0*x[r]*x[r]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2))-
                    2.0/((1.0+radius2)*(1.0+radius2));
         else
            value = 8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
      else
         if (r == s)
            value = -8.0*x[r]*x[r]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2))+
                     2.0/((1.0+radius2)*(1.0+radius2));
         else
            value = -8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
   
      return value;
   }
   else if (globtestnum == 1)
   {
      double value, radius2;
      double epsp = EPSILONP, epsm = EPSILONM;

      if (r == s)
         if (thesign < 0.0)
            value = 2.0*epsp/epsm;
         else
            value = 2.0;
      else
         value = 0.0;
   
      return value;
   }
   else if (globtestnum == 2)
   {
      double value, radius2;
      double epsp = EPSILONP, epsm = EPSILONM;
      int t;

      radius2 = 0.0;
      for (t = 0; t < grid.dim; t++)
         radius2 += x[t]*x[t];

      if (r == s)
         if (thesign < 0.0)
            value = (2.0*epsp/epsm+4.0*(epsp/epsm)*(epsp/epsm)*x[r]*x[r])*
                    exp(epsp/epsm*radius2-(epsp/epsm-1.0)*grid.radius*grid.radius);
         else
            value = (2.0+4.0*x[r]*x[r])*exp(radius2);
      else
         if (thesign < 0.0)
            value = 4.0*(epsp/epsm)*(epsp/epsm)*x[r]*x[s]*
                    exp(epsp/epsm*radius2-(epsp/epsm-1.0)*grid.radius*grid.radius);
         else
            value = 4.0*x[r]*x[s]*exp(radius2);
   
      return value;
   }
   else if (globtestnum == 3 || globtestnum == 4)
   {
      double value, radius2;
      double epsp = EPSILONP, epsm = EPSILONM;
      int t;
      double A = 1.0, B = fabs(1.0-epsp/epsm), 
             C = epsp/epsm, D = grid.radius*grid.radius*(1.0-epsp/epsm)+B;

      radius2 = 0.0;
      for (t = 0; t < grid.dim; t++)
         radius2 += x[t]*x[t];

      if (thesign < 0.0)
      {
         value = -C*x[r]*x[s]/(C*radius2+D);
         if (r == s)
            value += 1.0;
         value *= C/sqrt(C*radius2+D);
      }
      else
      {
         value = -A*x[r]*x[s]/(A*radius2+B);
         if (r == s)
            value += 1.0;
         value *= A/sqrt(A*radius2+B);
      }
   
      return value;
   }

   else if (globtestnum == 12 ||globtestnum == 11 || globtestnum == 10)
   {
      return getD2uTest(x,r,s,thesign);
   }

   cerr << "undefined getD2u" << endl;
   exit(-1);
   return 0.0;
}

double getD2u(int *index, int r, int s, int rstar, int sstar, double alpha, 
              double thesign, GridData &grid)
{
   double x[grid.dim];
   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];
   return getD2u(x, r, s, thesign,grid);
   
}


double geta(double *x, double thesign, PBData& pb, GridData& grid){
  if (globtesta == 1){
    // used in globtestnum = 11, 6 surface, does not change result
    assert(globtestnum == 11 && "globesta 1 goes with globtestnum 11");

    if (thesign < 0.0)
        return pb.epsilonm * sin(x[0]);
    else
        return pb.epsilonp * cos(x[grid.dim-1]);
  
  }else if (globtesta == 2){
    // used in globtestnum = 0, motion
    assert(globtestnum == 0 && "globesta 1 goes with globtestnum 0");
    double r2 = pow(x[0],2) + pow(x[1],2) + pow(x[2],2);
    if (thesign < 0.0)
      return pb.epsilonm * sin(r2);
    else
      return pb.epsilonp * cos(r2); 
  
  }else if(globtesta == 3){
    assert(globtestnum == 11 && "globesta 3 goes with globtestnum 11");
    if (thesign < 0.0)
        return pb.epsilonm * (sin(x[0]) + exp(x[1] * x[1]) + pow(x[2],3));
    else
        return pb.epsilonp * (cos(x[2]) + exp(- x[1] * x[1]) + pow(x[0],3));

  }else{
    return 0;
   }
   return 0;
}

void geta(double ***a, double ***S, PBData &pb, GridData &grid)
{
  for (int i = 0; i <= grid.nx[0]; i++){
    for (int j = 0; j <= grid.nx[1]; j++){
      for (int k = 0; k <= grid.nx[2]; k++){
        int tindex[3] = {i, j, k};
        double x[3];
        sub2coord(x,tindex,grid);
        double thesign = (S[i][j][k]<0.0)? -1.0: 1.0;
        a[i][j][k] = geta(x, thesign, pb, grid);
      }
    }
  }
}

double gettau(double *x, GridData &grid)
{
   if (globtestnum == 0)
   {
      int r;
      double valuep, valuem, radius2;

      radius2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
   
      valuep = -1.0/(1.0+radius2);
      valuem = 1.0/(1.0+radius2);
      return valuep-valuem;
   }
   else if (globtestnum == 1)
      return 0.0;
   else if (globtestnum == 2)
      return 0.0;
   else if (globtestnum == 3 || globtestnum == 4)
      return 0.0;
   else if (globtestnum == 9)
      return gettau_custom(x,grid);
   else if (globtestnum == 12 ||globtestnum == 11 || globtestnum == 10)
      return gettauTest(x);
   
      

   cout << "NOT SUPPOSED TO BE HERE" << endl;
   return 0.0;
}

void gettau(double &tau, int *index, int rstar, int sstar, double alpha, 
            GridData &grid)
{
   double x[grid.dim];
   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];
   tau = gettau(x, grid);
}

void getDtau(double *Dtau, double *x, GridData &grid)
{
   int r, s;
  
   if (globtestnum == 0)
   {
      double valuep, valuem, radius2;

      radius2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
   
      for (s = 0; s < grid.dim; s++)
      {
         valuep = 2.0*x[s]/((1.0+radius2)*(1.0+radius2));
         valuem = -2.0*x[s]/((1.0+radius2)*(1.0+radius2));
         Dtau[s] = valuep-valuem;
      }
   }
   else if (globtestnum == 1)
      for (s = 0; s < grid.dim; s++)
         Dtau[s] = 0.0;
   else if (globtestnum == 2)
      for (s = 0; s < grid.dim; s++)
         Dtau[s] = 0.0;
   else if (globtestnum == 3 || globtestnum == 4)
      for (s = 0; s < grid.dim; s++)
         Dtau[s] = 0.0;
   else if (globtestnum == 9)
      getDtau_custom(Dtau,x,grid);
    else if (globtestnum == 12 ||globtestnum == 11 || globtestnum == 10){
      getDtauTest(Dtau,x);
    }
}

void getDtau(double *Dtau, int *index, int rstar, int sstar, double alpha, 
             GridData &grid)
{
   double x[grid.dim];
   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];
   getDtau(Dtau, x, grid);
}

void getD2tau(double **D2tau, double *x, GridData &grid)
{
   int r, s, t;

   if (globtestnum == 0)
   {
      double valuep, valuem, radius2;

      radius2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
      for (r = 0; r < grid.dim; r++)
         for (s = r; s < grid.dim; s++)
         {
            valuep = -8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
            valuem = 8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
            if (s == r)
            {
               valuep += 2.0/((1.0+radius2)*(1.0+radius2));
               valuem += -2.0/((1.0+radius2)*(1.0+radius2));
            }
            D2tau[r][s] = valuep-valuem;
            D2tau[s][r] = D2tau[r][s];
         }
   }
   else if (globtestnum == 1)
   {
      for (r = 0; r < grid.dim; r++)
         for (s = r; s < grid.dim; s++)
         {
            D2tau[r][s] = 0.0;
            D2tau[s][r] = D2tau[r][s];
         }
   }
   else if (globtestnum == 2){
      for (r = 0; r < grid.dim; r++)
         for (s = r; s < grid.dim; s++)
         {
            D2tau[r][s] = 0.0;
            D2tau[s][r] = D2tau[r][s];
         }
   }
   else if (globtestnum == 3 || globtestnum == 4){
      for (r = 0; r < grid.dim; r++)
         for (s = r; s < grid.dim; s++)
         {
            D2tau[r][s] = 0.0;
            D2tau[s][r] = D2tau[r][s];
         }
   }
   else if (globtestnum == 9){
    getD2tau_custom(D2tau, x, grid);
  }
  else if (globtestnum == 12 ||globtestnum == 11 || globtestnum == 10){
    getD2tauTest(D2tau, x);
  }
  else{
   cerr << "undefiend sigma" << endl;
   exit(-1); 
  }
  
}

void getD2tau(double **D2tau, int *index, int rstar, int sstar, double alpha, 
              GridData &grid)
{
   double x[grid.dim];
   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];
   getD2tau(D2tau, x, grid);

}

double getsigma(double *x, double *normal, PBData &pb, GridData &grid)
{
   if (globtestnum == 0)
   {
      int r, s;
      double sigma, gradp[grid.dim], gradm[grid.dim], radius2;
   
      sigma = 0.0;
      radius2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
      for (r = 0; r < grid.dim; r++)
      {
         gradp[r] = 2.0*x[r]/((1.0+radius2)*(1.0+radius2));
         gradm[r] = -2.0*x[r]/((1.0+radius2)*(1.0+radius2));
         sigma += (pb.epsilonp*gradp[r]-pb.epsilonm*gradm[r])*normal[r];
      }

      return sigma;
   }
   else if (globtestnum == 1)
      return 0.0;
   else if (globtestnum == 2)
      return 0.0;
   else if (globtestnum == 3 || globtestnum == 4)
      return 0.0;
   else if (globtestnum == 9){
      return getsigma_custom(x,normal,pb,grid);
    }
    else if (globtestnum == 12 ||globtestnum == 11 || globtestnum == 10){
      return getsigmaTest(x,normal,pb);
    }

   cerr << "undefiend sigma" << endl;
   exit(-1);
   return 0.0;
}

void getsigma(double &sigma, int *index, int rstar, int sstar, double alpha, 
              double *normal, PBData &pb, GridData &grid)
{
   double x[grid.dim];
   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];
   sigma = getsigma(x,normal,pb,grid);
}

void getDsigma(double *Dsigma, double *x, double *normal, double **Dnormal, PBData &pb, 
               GridData &grid)
{
   int r, s, t;
  
   if (globtestnum == 0)
   {
      double gradp[grid.dim], gradm[grid.dim], D2p[grid.dim][grid.dim],
             D2m[grid.dim][grid.dim], radius2;

      radius2 = 0.0;
      for (r = 0; r < grid.dim; r++)
         radius2 += x[r]*x[r];
      for (r = 0; r < grid.dim; r++)
      {
         gradp[r] = 2.0*x[r]/((1.0+radius2)*(1.0+radius2));//grad(u+)
         gradm[r] = -2.0*x[r]/((1.0+radius2)*(1.0+radius2));//grad(u-)
         for (s = 0; s < grid.dim; s++)
         {
            D2p[r][s] = -8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
            D2m[r][s] = 8.0*x[r]*x[s]/((1.0+radius2)*(1.0+radius2)*(1.0+radius2));
            if (s == r)
            {
               D2p[r][s] += 2.0/((1.0+radius2)*(1.0+radius2));
               D2m[r][s] += -2.0/((1.0+radius2)*(1.0+radius2));
            }
         }
      }
      for (r = 0; r < grid.dim; r++)
      {
         Dsigma[r] = 0.0;
         for (s = 0; s < grid.dim; s++)
            Dsigma[r] += (pb.epsilonp*gradp[s]-pb.epsilonm*gradm[s])*Dnormal[s][r]+
                         (pb.epsilonp*D2p[r][s]-pb.epsilonm*D2m[r][s])*normal[s];
      }
   }
   else if (globtestnum == 1)
      for (r = 0; r < grid.dim; r++)
         Dsigma[r] = 0.0;
   else if (globtestnum == 2)
      for (r = 0; r < grid.dim; r++)
         Dsigma[r] = 0.0;
   else if (globtestnum == 3 || globtestnum == 4)
      for (r = 0; r < grid.dim; r++)
         Dsigma[r] = 0.0;
   else if (globtestnum == 9){
      getDsigma_custom(Dsigma, x, normal, Dnormal, pb, grid);
    }
    else if (globtestnum == 12 ||globtestnum == 11 || globtestnum == 10){
      getDsigmaTest(Dsigma, x, normal, Dnormal, pb);
    }
}

void getDsigma(double *Dsigma, int *index, int rstar, int sstar, double alpha, 
              double *normal, double **Dnormal, PBData &pb, GridData &grid)
{

      double x[grid.dim];
      sub2coord(x,index,grid);
      x[rstar] += sstar*alpha*grid.dx[rstar];
      getDsigma(Dsigma, x, normal, Dnormal, pb, grid);
}



void init_surf(double ***S, double radius, GridData &grid, int opt){
   double x[3];
   
   if(opt == 0){
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               int index[3] = {i,j,k};
               sub2coord(x,index,grid);
               S[i][j][k] = sqrt(pow(x[0],2) + pow(x[1],2) + pow(x[2],2)) - radius;
            }
         }
      }
   }else if (opt == 1){
      // strange torus from getinit
      double radii[3];
      radii[0] = 0.5;
      radii[1] = 0.25;
      radii[2] = 0.15;
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               int index[3] = {i,j,k};
               sub2coord(x,index,grid);

               double value = 0;
               for (int r = 1; r < grid.dim; r++)
                  value += x[r]*x[r];
               value = sqrt(value)-radii[0];
               double theta = atan2(x[2],x[1]);
   // main one used in tests
               value = sqrt(value*value+x[0]*x[0])-(0.5*(radii[1]-radii[2])*sin(3.0*theta)+
                                              0.5*(radii[1]+radii[2]));
   //      value = sqrt(value*value+x[0]*x[0])-(0.5*(radii[1]-radii[2])*sin(theta)+
   //                                           0.5*(radii[1]+radii[2]));
   //      value = sqrt(value*value+x[0]*x[0])-radii[1];
               S[i][j][k] = value;
            }
         }
      }

   }
   else if (opt == 11){
      //ellipsoid
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               int index[3] = {i,j,k};
               sub2coord(x,index,grid);
               S[i][j][k]= 2 * pow(x[0],2) + 3 * pow(x[1],2) + 6 * pow(x[2],2) - 1.3 * 1.3;
            }
         }
      }
   }else if (opt == 12){
      //donut from paper
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               int index[3] = {i,j,k};
               sub2coord(x,index,grid);
               S[i][j][k]= -0.3 * 0.3 + pow(-0.6 + sqrt(pow(x[0],2) + pow(x[1],2)),2) + pow(x[2],2);
            }
         }
      }
   }else if (opt == 13){
      //banana
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               int index[3] = {i,j,k};
               sub2coord(x,index,grid);
               S[i][j][k] = 1521 - 94*pow(6 + 7*x[0],2) + pow(6 + 7*x[0],4) + 3822*pow(x[1],2) + 2401*pow(x[1],4) - 4606*pow(x[2],2) + 4802*pow(x[1],2)*pow(x[2],2) + 3601.5*pow(x[2],4) + 98*pow(6 + 7*x[0],2)*(pow(x[1],2) + pow(x[2],2));
            }
         }
      }
   }else if (opt == 14){
      //8 balls
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               int index[3] = {i,j,k};
               sub2coord(x,index,grid);
               double d = std::numeric_limits<double>::infinity();
               for(int n = 0; n <=7; n++){
                  double c[3] = {pow(-1,floor(n/4)) * 0.5, pow(-1,floor(n/2)) * 0.5, pow(-1,n) * 0.5};
                  d = min(d,dist(x,c)-0.3);
                  S[i][j][k] = d;
               }
            }
         }
      }
   }else if (opt == 15){
      //peanut
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               int index[3] = {i,j,k};
               sub2coord(x,index,grid);
               double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

               double theta = atan2(x[1],x[0]);
               
               double phi = acos(x[2]/r);

               if(r < 1e-9){
                S[i][j][k] = 0 - 0.5 - 0.2;
               }else{
                S[i][j][k] = r - 0.5 - 0.2 * sin(2*theta)*sin(phi); 
               }               
               
            }
         }
      }
   }
   else if (opt == 16){
      //popcorn
      double r0 = 0.6;
      double A = 2;
      double sigma = 0.2;
      double xk[12][3] ;

      for(int k = 0; k <=11; k++){
         if(k<=9){
            xk[k][0] = r0/sqrt(5.0) * 2 * cos(2 * k * M_PI/5 - floor(k/5.0) * M_PI);
            xk[k][1] = r0/sqrt(5.0) * 2 * sin(2 * k * M_PI/5 - floor(k/5.0) * M_PI);
            xk[k][2] = r0/sqrt(5.0) * pow(-1,floor(k/5.0));
         }else{
            xk[k][0] = r0 * 0;
            xk[k][1] = r0 * 0;
            xk[k][2] = r0 * pow(-1,k-10);
         }
      }
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               int index[3] = {i,j,k};
               sub2coord(x,index,grid);
               S[i][j][k] = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) - r0;
               double temp = 0.0;
               for(int n = 0; n <=11; n++){
                  temp += 2 * exp( -25 * pow(dist(x,xk[n]),2));
               }
               S[i][j][k] -= temp;
            }
         }
      }
   }else if (opt == 20){
      init_surf_protein_paper(S, grid);
   }
   else{
      cerr<<"invalid surface option";
      exit(1);
   }
}





void read_protein(double ***S, GridData &grid, vector< vector<double> >& p){
  cout<<"using protein file "<<PROTEINFILE<<endl;
   ifstream infile( PROJECTDIR + PROTEINFILE);
   if(!infile.is_open()){
      cerr<<"file do not exist"<<endl;
      exit(1);
   }

   string line;
   while(getline(infile,line)){
      stringstream ss(line);
      string x,y,z,r;
      getline(ss,x,',');
      getline(ss,y,',');
      getline(ss,z,',');
      getline(ss,r,',');
   
      
      p.push_back(std::vector<double> {stod(x),stod(y),stod(z),stod(r)});
   }

   for(int i = 0; i < 3; i++){
      printf("%d %f %f %f %f\n",i,p[i][0],p[i][1],p[i][2],p[i][3]);
   }

   for(int i = p.size()-3; i < p.size(); i++){
      printf("%d %f %f %f %f\n",i,p[i][0],p[i][1],p[i][2],p[i][3]);
   }   
}


inline double chi(double x){
   const double eta = 1.0/ETA; //in paper 1/40
   return 1.0/2.0 * (1.0 + tanh(x/eta));
}

void init_surf_protein_paper(double ***S, GridData& grid){
  // if file exist
  string file = PROTEINDIR + string("surf_protein_paper") + "_eta"+to_string( (int)ETA) + "_gd" + to_string(GRIDNUM) + "_" + PROTEINFILE;
  if(fileExist(file)){
    cout<<"reading "<< file<<endl;
    readMatrix(file, S, grid);
    return;
  }
  cout<<"file "<<file <<" do not exist, construct from protein"<<endl;
  // if file do not exist
  vector< vector<double> > p;
  read_protein(S,grid,p);
  const double c = 0.25;
   #pragma omp parallel for collapse(3)
   for(int i = 0; i <= grid.nx[0]; i++){
      for(int j = 0; j <= grid.nx[1]; j++){
         for(int k = 0; k <= grid.nx[2]; k++){
            array<double,3> x = sub2coord(array<int,3> {i,j,k}, grid);
            double sum = 0.0;

            for(int n = 0; n < p.size(); n++){
               array<double,3> px = {p[n][0],p[n][1],p[n][2]};//px is coordinate
               sum += chi(p[n][3] - dist(x,px)); //p[n][3] = radius of atom n
            }
            S[i][j][k] = c - sum;
         }
      }
   }

   // FillVoid(S, grid);
   // write to file
   cout<<"write "<< file<<endl;
   write_field(file, S, grid);

}



// get maximum error of u at grid points
void checkanswer(double ***u, double ***S,GridData &grid)
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

      // if (globwriteerr){
      //   outfile_uerr<<tindex[0]<<","<<tindex[1]<<","<<tindex[2]<<",";
      //   outfile_uerr <<setprecision(12)<<tmperr<< endl;;
      // }


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



// if surface on grid point, within tol, shift to tol
void perturbelt(double ***S, int *index, double tol)
{
   if (evalarray(S,index) >= 0.0 && evalarray(S,index) < tol)
      setvalarray(S,index,tol);
   else if (evalarray(S,index) < 0.0 && evalarray(S,index) > -tol)
      setvalarray(S,index,-tol);
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






// compute exact radius at thetime for motion test, radius0 is initial radius,
// x0 is radius at current step
// some are computed exactly, some solve ODE by one step of RK
double getexactradius(double thetime, double radius0, double x0, double tol, int Nstep,
                      GridData &grid)
{
   // find root directly using newton's method
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
   // analytic solution of r(t)
   else if (globtestnum == 1)
   {
      double epsp = EPSILONP, epsm = EPSILONM;
      return grid.radius0*exp(2.0*(1.0-epsp/epsm)*thetime);
   }
   // one step RK
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
   // one step RK from 
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
                double **jumpuxxcoef,
               int *index, int rstar, int sstar, double alpha, 
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




double evalcoef(double u0, double ***ucoef, double *uxcoef,
                double *uxxcoef, 
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



//uxx include cross derivative
double evalcoef(double u0, double ***ucoef, double *uxcoef, 
                double **uxxcoef, 
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


// get the exact maximum vn for all time, used in cfl condition
double getexactvnmax(PBData &pb){
   double exactmaxvn;
   if (globtestnum == 0)
      exactmaxvn = 3.0*sqrt(3.0)/4.0;
   else if (globtestnum == 1)
      exactmaxvn = fabs(2.0*(1.0-pb.epsilonp/pb.epsilonm)); //taken at r = 1, 429.48
   else if (globtestnum == 2)
      exactmaxvn = 2.0*exp(1.0)*fabs(1.0-pb.epsilonp/pb.epsilonm); //taken at r = 1, 158
   else if (globtestnum == 3 || globtestnum == 4)
      exactmaxvn = fabs(1.0-pb.epsilonp/pb.epsilonm)/
                   sqrt(1.0+fabs(1.0-pb.epsilonp/pb.epsilonm));
   return exactmaxvn;
}


void setgriddt(PBData &pb, GridData &grid){
   double maxvn = getexactvnmax(pb);
   
   if (!globdtdx2)
      grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim*pb.gamma0+2.0*maxvn*grid.mindx);
   else if (globtestnum == 0);
   else if (globtestnum == 1)
      grid.dt = grid.mindx*grid.mindx/22.0;
   else if (globtestnum == 2)
      grid.dt = grid.mindx*grid.mindx/58.0;
   else if (globtestnum == 3 || globtestnum == 4)
      grid.dt = grid.mindx*grid.mindx/1.2;
   else{
      printf("undefined globtestnum");
      exit (0);
   }

   printf("dt = %f\n", grid.dt);
}

void checkwithexactvn(double ***vn, double ***S, PBData &pb, GridData &grid){
   double exactvn;
   if (globtestnum == 0){
      exactvn = 4.0*grid.radius/((1.0+grid.radius*grid.radius)*
                                 (1.0+grid.radius*grid.radius));
   }
   else if (globtestnum == 1){
      exactvn = 2.0*(1.0-pb.epsilonp/pb.epsilonm)*0.5*
                exp(2.0*(1.0-pb.epsilonp/pb.epsilonm)*grid.t);
   }
   else if (globtestnum == 2){
      exactvn = 2.0*grid.radius*exp(grid.radius*grid.radius)*
                (1.0-pb.epsilonp/pb.epsilonm);
   }
   else if (globtestnum == 3 || globtestnum == 4){
      exactvn = (1.0-pb.epsilonp/pb.epsilonm)*grid.radius/
                sqrt(grid.radius*grid.radius+fabs(1.0-pb.epsilonp/pb.epsilonm));
   }
   else{
      cerr<<"invalid motion test"<<endl;
      exit(1);
   }      


   vector<double> vnAll;
   vector<double> vnInterface;
   array<double,3> xAll; // coordinate of max vn extended
   array<double,3> xInterface; // coordinate of max vn interface
   double maxvnAll = 0; // max vn extended
   double maxvnInterface = 0; // max vn interface
   for(int i = 0; i <= grid.nx[0]; i++){
      for(int j = 0; j <= grid.nx[1]; j++){
         for(int k = 0; k <= grid.nx[2]; k++){
            array<int,3> x = {i, j, k};
            double absvn = abs(vn[i][j][k]);
            vnAll.push_back(absvn);
            if(absvn > maxvnAll){
               maxvnAll = absvn;
               xAll = sub2coord(array<int,3>{i,j,k},grid);
            }
            
            if(nearinterface(x, S, grid)){
               vnInterface.push_back(absvn);
               if(absvn > maxvnInterface){
                  maxvnInterface = absvn;
                  xInterface = sub2coord(array<int,3>{i,j,k},grid);
               }

            }
         }
      }
   }
   printf("abs(exactvn) = %f\n", abs(exactvn));
   printf("vn interface mean %f, variance %f, rmse %f\n", mean(vnAll), variance(vnAll), rmse(vnAll, abs(exactvn)));
   printf("maxvn interface %f at (%f,%f,%f)\n",maxvnInterface, xInterface[0], xInterface[1], xInterface[2]);
   printf("vn extend mean %f, variance %f, rmse %f\n", mean(vnInterface), variance(vnInterface), rmse(vnInterface, abs(exactvn)));
   printf("maxvn extend %f at (%f,%f,%f)\n",maxvnAll, xAll[0], xAll[1], xAll[2]);
   
}