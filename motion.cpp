#include "helper.h"
#include "ccim.h"
#include "cim12.h"
#include "icim.h"
#include "input.h"
#include "solver.h"
#include "hypresolver.h"
#include "march.h"
#include "advance.h"

using namespace std;


void solvepde(double*** S, double ***a, PBData &pb, GridData &grid){
	TempStruct tmp;//sparse matrix solver
	init_TempStruct(tmp, grid);
	
	StorageStruct *Dusmall;
	int smallsize = 0;


	linearsystemcim345(tmp.A,tmp.b,Dusmall,smallsize,a,S,pb,grid);

	HypreSolve(pb.psi, tmp.A, tmp.b, grid, S, pb);

	clearfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
}


void rhsvngrad(double*** rhs, double*** u, MarchStruct &march, PBData &pb, GridData &grid){
	int r, i, s;
	int tindex[3], rindex[3];
	int deg = 3;
	double u1d[2*deg+1];
	double dfp[grid.dim], dfn[grid.dim], df[grid.dim];
	double grad;
	double maxvn = 0;

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
}


void eulerstep(double*** u, double ***rhs, double dt, GridData &grid){
	int i;
	int tindex[grid.dim];
	for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(u,tindex,evalarray(u,tindex)+grid.dt*evalarray(rhs,tindex));
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
}

int main(int argc, char* argv[])
{

	cmdLine(argc, argv);
	
	double tol = 1.0e-12; // tol for pertubation. tol for solver in grid

	GridData grid;
	init_grid(grid, GRIDNUM, 1.0);

	
	PBData pb;
	init_PBData(pb, grid, EPSILONP, EPSILONM);

	// create initial surface
	double ***S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	init_surf(S, grid.radius, grid, SURFOPT);

	double ***a = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

	double ***rhs = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	
	StorageStruct *Dusmall;
	int smallsize = 0;

	MarchStruct march;
	init_march(march, S, pb, grid);

	march.Dusmall = Dusmall;
	march.smallsize = smallsize;
	march.dorig = S;


	for (int step = 1; step <= runStep && grid.tfinal-grid.t > 1.0e-14; step++)
	{
		
		cout << "\nSTEP = " << step << " and time = " << grid.t << endl;
		
		if (globperturb == 0)
			perturb(S,grid.tol,pb,grid);

		solvepde(S, a, pb, grid);

		fmarch(march,-1.0,grid);

		rhsvngrad(rhs, S, march, pb, grid);

		eulerstep(S, rhs, grid.dt, grid);


		cout << "Exact radius = " << grid.radius << endl;
		checkwithexact(S, grid.radius, grid); // check radius with exact
		
	}	

	

	
	free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
	free_matrix(a,grid.nx[0],grid.nx[1],grid.nx[2]);
	free_matrix(rhs,grid.nx[0],grid.nx[1],grid.nx[2]);
}




