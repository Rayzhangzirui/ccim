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
