#include "tryvn.h"
#include "helper.h"
#include "advance.h"
#include "global.h"

using namespace std;


int main(int argc, char* argv[])
{

	cmdLine(argc, argv);
	
	string file_prefix = string("data") + "_n"  + to_string(GRIDNUM)
										 + "_c"  + to_string(globcim)
										 + "_s"  + to_string(SURFOPT)
										 + "_t"  + to_string(globtestnum)
										 + "_a"  + to_string(globtesta)
										 + "_d"  + to_string(globdist)
										 + "_epsm"  + to_string((int)round(EPSILONM));
										 
	globwritemx = false;
	globwriteerr = false;								 
	if(globwritemx){
		outfile_mx.open(file_prefix + "_mx.dat",ofstream::out);	
	}
	if(globwriteerr){
		
		outfile_uerr.open(file_prefix + "_uerr.dat",ofstream::out);
		outfile_Duerr.open(file_prefix + "_Duerr.dat",ofstream::out);
	}
	
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
	else if (globperturb > 0)
		perturbstatus(S,grid.tol,globperturb,pb,grid);

	MarchStruct march;
	init_march(march, S, pb, grid);

	
	// if (globtestnum == 0)
	// 	grid.tfinal = 15.0*2.0/50.0/(2.0*3.0*sqrt(3.0)/4.0);
	// else if (globtestnum == 1)
	// 	grid.tfinal = 10.0*2.0/50.0/(2.0*2.0*fabs(1.0-pb.epsilonp/pb.epsilonm));
	// else if (globtestnum == 2)
	// 	grid.tfinal = 10.0*2.0/50.0/(2.0*2.0*exp(1.0)*fabs(1.0-pb.epsilonp/pb.epsilonm));
	// else if (globtestnum == 3 || globtestnum == 4){
	// // grid.tfinal = 10.0*2.0/50.0/(2.0*fabs(1.0-pb.epsilonp/pb.epsilonm)/
	// // 											sqrt(1.0+fabs(1.0-pb.epsilonp/pb.epsilonm)));
	// 		grid.tfinal = 0.1;
	// }

	grid.tfinal = globtime;

	double ***a = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

	// check initial error in radius
	checkwithexact(S, grid.radius,grid);
	for (int step = 1; step <= runStep && grid.tfinal-grid.t > 1.0e-14; step++)
	{
		
		cout << "\nSTEP = " << step << " and time = " << grid.t << endl;
		grid.maxvn = 0.0;
		if (globtesta == 0)
			advance(S,pb,march,tmp,grid);
		else
		{
			geta(a,pb.psi,S,pb,grid);
			advance(S,a,pb,march,tmp,grid);
		}
		
		if (globperturb == 0)
			perturb(S,grid.tol,pb,grid);
		else if (globperturb > 0)
			perturbstatus(S,grid.tol,globperturb,pb,grid);
		cout << "Exact radius = " << grid.radius << endl;
		checkwithexact(S, grid.radius, grid); // check radius with exact
		
	}
	cout << "Final time = " << grid.t << endl;
	clearfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);

	write_field( "final_motion_surf_t" + to_string(globtestnum) + "_gn" + to_string(GRIDNUM) + ".txt", S, grid);//output final surface

	free_matrix(a,grid.nx[0],grid.nx[1],grid.nx[2]);
	free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
}




