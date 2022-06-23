#include "pb.h"

void init_PBData(PBData &pb, GridData &grid, double epsp, double epsm)
{
   pb.epsilonm = epsm;
   pb.epsilonp = epsp;
   pb.gamma0 = 0.0;
   pb.psi = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   pb.dim = grid.dim;
   pb.N = 2;//no solute
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
   for (int i = 0; i < pb.N; i++)
   {
      pb.Q[i] = 0.0;//1.0;
      pb.c[i] = 0.0;//1.0;
      pb.epsilone[i] = 0.0;//1.0;
      pb.epsilonlj[i] = 0.0;//0.0159;
      pb.sigmalj[i] = 0.0;//3.653;
   }
   pb.beta = 1.0;
   pb.rho0 = 0.0333;
   pb.LJ = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);


   for(int i = 0; i <= grid.nx[0]; i++){
      for(int j = 0; j <= grid.nx[1]; j++){
         for(int k = 0; k <= grid.nx[2]; k++){
            pb.psi[i][j][k] = 0;
            pb.LJ[i][j][k] = 0;
         }
      }
   }

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