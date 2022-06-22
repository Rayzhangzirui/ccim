#include "helper.h"
#include "ccim.h"
#include "cim12.h"
#include "icim.h"
#include "input.h"
#include "solver.h"
#include "hypresolver.h"

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
		perturb(S, grid.tol, pb, grid);
	// else if (globperturb > 0)
	// 	perturbstatus(S,grid.tol,globperturb,pb,grid);

	double ***a = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

	StorageStruct *Dusmall;
	int smallsize = 0;
	
	if (globcim == 7){
		double ***sign_surf = matrix(grid.nx[0],grid.nx[1],grid.nx[2]); // sign of flipped surface
		double *** psi_true = matrix(grid.nx[0],grid.nx[1],grid.nx[2]); // reconstructed psi

		SignSurf(sign_surf, S, grid);
		int maxIter = 10;
		int depth = 1;
		bool flipbydet = true;
		flip(sign_surf, S, pb, grid, maxIter, depth, flipbydet);
		linearsystem_icim(tmp.A,tmp.b, S, sign_surf, pb, grid);


		HypreSolve(pb.psi, tmp.A, tmp.b, grid, S, pb);

		// reconstruct solution
		reconstruct(psi_true, S, pb.psi, sign_surf, grid);

		CheckErrGrid(psi_true, S, pb, grid); // check original surface 
		CheckIcimDu(pb.psi, psi_true, sign_surf, S, pb, grid);

		free_matrix(sign_surf, grid.nx[0],grid.nx[1],grid.nx[2]);
		free_matrix(psi_true, grid.nx[0],grid.nx[1],grid.nx[2]);

    }
    if (globcim == 345) {
		linearsystemcim345(tmp.A,tmp.b,Dusmall,smallsize,a,S,pb,grid);

		HypreSolve(pb.psi, tmp.A, tmp.b, grid, S, pb);
		checkanswer(pb.psi,S,grid);
		checkDuStorage(pb.psi, Dusmall, smallsize, S, pb, grid);
    }
	

	

	clearfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
	free_matrix(a,grid.nx[0],grid.nx[1],grid.nx[2]);
	free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
}




