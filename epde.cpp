#include "helper.h"
#include "ccim.h"
#include "input.h"
#include "solver.h"

using namespace std;

int main(int argc, char* argv[])
{

	cmdLine(argc, argv);
	
	double tol = 1.0e-14; // tol for pertubation. tol for solver in grid
	double gridbound = 1.0;

	GridData grid;
	init_grid(grid, GRIDNUM, gridbound);

	TempStruct tmp;//sparse matrix solver
	init_TempStruct(tmp, grid);
	
	PBData pb;
	init_PBData(pb, grid, EPSILONP, EPSILONM);

	// create initial surface
	double ***S;
	S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	init_surf(S, grid.radius, grid, SURFOPT);

	if (globperturb == 0)
		perturb(S,grid.tol,pb,grid);
	// else if (globperturb > 0)
	// 	perturbstatus(S,grid.tol,globperturb,pb,grid);

	double ***a = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

	StorageStruct *Dusmall;
	int smallsize = 0;

	linearsystemcim345(tmp.A,tmp.b,Dusmall,smallsize,a,S,pb,grid);

  	BICGSTABsmall(pb.psi,tmp.A,tmp.b,grid,grid.nx[0]*grid.nx[1]*grid.nx[2],tollinsolve,S,pb,tmp);
	// amgsolve(pb.psi, tmp.A, tmp.b, grid, S, pb);

	checkanswer(pb.psi,S,grid);
	checkDuStorage(pb.psi, Dusmall, smallsize, S, pb, grid);

	clearfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);

	free_matrix(a,grid.nx[0],grid.nx[1],grid.nx[2]);
	free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
}




