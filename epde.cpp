#include "helper.h"
#include "ccim.h"
#include "cim12.h"
#include "icim.h"
#include "input.h"
#include "solver.h"
#include "amg.h"
#include "hypresolver.h"

using namespace std;


void LinearSolver(double ***x, SparseElt2**** &A, double ***b, double ***a, double ***u,  GridData &grid, 
                   int numsteps, double tol, PBData &pb, 
                   TempStruct &tmp)
{
    clock_t cstart = clock ();

    char name[20];

    switch (globlinsolve){
        case 0:
            strcpy(name, "bicgstab");
            BICGSTABsmall(x,A, b, grid, numsteps, tollinsolve, a, u, pb,tmp);
            break;
        case 1:
            strcpy(name, "ilu-bicgstab");
            ILUsmall(tmp.M, A, grid, a, u,pb);
      
            prerightBICGSTABsmall(x, A, b, tmp.M, grid, numsteps, tollinsolve, a, u, pb,tmp);
          
            gaussseidelsmall(x,A, b, grid,globGSsmooth,a,u,pb);
            break;
        case 2:
            strcpy(name, "my-amg");
            AMGsmall3(x, A, b, grid, 2, 0, grid.nx[0], numsteps, tollinsolve, a, u, pb);
            break;
        case 4:
            strcpy(name, "hypre-amg");
            HypreSolve(x, A, b,  grid, u, pb);
            break;
        case 5:
            strcpy(name, "hypre-ilu-bicgstab");
            HypreSolve(x, A, b,  grid, u, pb);
            break;


        default:
            cerr << "No linear solve option" << endl;
            exit(1);
    }
    clock_t cend = clock ();
    printf("%s clock time = %f (s).\n", name, (double) (cend-cstart)/CLOCKS_PER_SEC );
}

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
    //  perturbstatus(S,grid.tol,globperturb,pb,grid);

    double ***a = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

    StorageStruct *Dusmall;
    int smallsize = 0;
    
    if (globcim == 7){
        double *** sign_surf = matrix(grid.nx[0],grid.nx[1],grid.nx[2]); // sign of flipped surface
        double *** psi_true = matrix(grid.nx[0],grid.nx[1],grid.nx[2]); // reconstructed psi

        SignSurf(sign_surf, S, grid);
        int maxIter = 10;
        int depth = 1;
        bool flipbydet = true;
        flip(sign_surf, S, pb, grid, maxIter, depth, flipbydet);
        linearsystem_icim(tmp.A,tmp.b, S, sign_surf, pb, grid);


        LinearSolver(pb.psi, tmp.A, tmp.b, a, S, grid, grid.N, tollinsolve, pb,tmp);

        // reconstruct solution
        reconstruct(psi_true, S, pb.psi, sign_surf, grid);

        CheckErrGrid(psi_true, S, pb, grid); // check original surface 
        CheckIcimDu(pb.psi, psi_true, sign_surf, S, pb, grid);

        free_matrix(sign_surf, grid.nx[0],grid.nx[1],grid.nx[2]);
        free_matrix(psi_true, grid.nx[0],grid.nx[1],grid.nx[2]);

    }

    if (globcim == 345 || globcim == 6) {
        linearsystemcim345(tmp.A,tmp.b,Dusmall,smallsize,a,S,pb,grid);

        LinearSolver(pb.psi, tmp.A, tmp.b, a, S, grid, grid.N, tollinsolve, pb,tmp);
        
        checkanswer(pb.psi,S,grid);
        checkDuStorage(pb.psi, Dusmall, smallsize, S, pb, grid);
    }



    
    

    clearfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
    free_matrix(a,grid.nx[0],grid.nx[1],grid.nx[2]);
    free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
}




