#include "sparse.h"
#include "helper.h"
#include "global.h"
#include "hypresolver.h"

int printLevel = 0;
/*
AMG solver, copy from ex5
*/
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
using namespace std;


void Copy3DArrayToHypre(double ***sx, HYPRE_IJVector &x, GridData &grid){
   for(int i = 0; i <= grid.nx[0]; i ++){
      for(int j = 0; j <= grid.nx[1]; j ++){
         for(int k = 0; k <= grid.nx[2]; k ++){
            Index ind = {i,j,k};
            int row = sub2ind(ind.data(), grid.nx, grid.dim); //row number in matrix
            
            HYPRE_IJVectorSetValues(x, 1, &row, &sx[i][j][k]);
            
         }
      }
   }
}


void Copy3DArrayFromHypre(double ***sx, HYPRE_IJVector &x, GridData &grid){
   for(int i = 0; i <= grid.nx[0]; i ++){
      for(int j = 0; j <= grid.nx[1]; j ++){
         for(int k = 0; k <= grid.nx[2]; k ++){
            Index ind = {i,j,k};
            int row = sub2ind(ind.data(), grid.nx, grid.dim); //row number in matrix
            
            HYPRE_IJVectorGetValues(x, 1, &row, &sx[i][j][k]);
         }
      }
   }
}




void HypreSolve(double ***sx, SparseElt2**** &sA, double ***sb, GridData &grid, double ***S, PBData &pb)
{
   int i;
   int myid, num_procs;
   int N, n;

   int ilower, iupper;
   int local_size, extra;

   
   int vis, print_system;


   HYPRE_IJMatrix A;
   HYPRE_ParCSRMatrix parcsr_A;
   HYPRE_IJVector b;
   HYPRE_ParVector par_b;
   HYPRE_IJVector x;
   HYPRE_ParVector par_x;
   

   HYPRE_Solver solver, precond;

   /* Initialize MPI */
   MPI_Init(NULL,NULL);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   /* Default problem parameters */
   
   
   vis = 0;
   print_system = 0;


   /* Preliminaries: want at least one processor per row */
   N = (grid.nx[0]+1) * (grid.nx[1]+1) * (grid.nx[2]+1); /* global number of rows */
   
   if (N < num_procs) n = sqrt(num_procs) + 1;
   

   /* Each processor knows only of its own rows - the range is denoted by ilower
      and upper.  Here we partition the rows. We account for the fact that
      N may not divide evenly by the number of processors. */
   local_size = N/num_procs;
   extra = N - local_size*num_procs;

   ilower = local_size*myid;
   ilower += hypre_min(myid, extra);

   iupper = local_size*(myid+1);
   iupper += hypre_min(myid+1, extra);
   iupper = iupper - 1;

   /* How many rows do I have? */
   local_size = iupper - ilower + 1;

   /* Create the matrix.
      Note that this is a square matrix, so we indicate the row partition
      size twice (since number of rows = number of cols) */
   HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

   /* Choose a parallel csr format storage (see the User's Manual) */
   HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

   /* Initialize before setting coefficients */
   HYPRE_IJMatrixInitialize(A);

    /* Create the rhs and solution */
   HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&b);
   HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(b);

   HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&x);
   HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
   HYPRE_IJVectorInitialize(x);


   // set up A, copy from sparse A
   {
      for(int i = 0; i <= grid.nx[0]; i ++){
         for(int j = 0; j <= grid.nx[1]; j ++){
            for(int k = 0; k <= grid.nx[2]; k ++){
               Index index = {i,j,k};
               
               double ehere;
               if (evalarray(S,index) < 0.0)
                  ehere = pb.epsilonm;
               else
                  ehere = pb.epsilonp;

               int row = sub2ind(index.data(), grid.nx, grid.dim); //row number in matrix
               
               HYPRE_IJVectorSetValues(x, 1, &row, &sx[i][j][k]); // at interior points, initial x is 0
               HYPRE_IJVectorSetValues(b, 1, &row, &sb[i][j][k]); // 
               
               

               if(atbound(index,grid)){                  
                  // set value at boundary
                  double value = 1;
                  
                  int ncols = 1;
                  HYPRE_IJMatrixAddToValues(A, 1, &ncols, &row, &row, &value);

                  HYPRE_IJVectorSetValues(x, 1, &row, &sb[i][j][k]); // set exact


               }else if (evalarray(sA,index.data()) == NULL){
                  // interior points
                  double value = 0;
                  for (int m = 0; m < grid.dim; m++){
                     value += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
                  }
                  
                  int ncols = 1;
                  // set diagonal
                  HYPRE_IJMatrixAddToValues(A, 1, &ncols, &row, &row, &value);
                  
                  // set off diagonal
                  for (int m = 0; m < grid.dim; m++){
                     for(int n : {-1,1}){
                        Index rindex = index;
                        rindex[m] = index[m] + n;
                        
                        int col = sub2ind(rindex.data(),grid.nx,grid.dim);
                        int ncols = 1;
                        value = -ehere/(grid.dx[m]*grid.dx[m]);
                        HYPRE_IJMatrixAddToValues(A, 1, &ncols, &row, &col, &value);
                     }
                  }

                  // set x and b
                  
               
               }else{
                  // interface points
                  SparseElt2 *current2;
                  for (current2 = evalarray(sA,index.data()); current2 != NULL; current2 = (*current2).next){
                     int col = sub2ind((*current2).cindex,grid.nx,grid.dim);
                     int ncols = 1;
                     
                     HYPRE_IJMatrixAddToValues(A, 1, &ncols, &row, &col, &((*current2).val));
                  }
               }


            }
         }
      }


   }

   

   /* Assemble after setting the coefficients */
   HYPRE_IJMatrixAssemble(A);

   /* Get the parcsr matrix object to use */
   HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);

   HYPRE_IJVectorAssemble(b);
   HYPRE_IJVectorGetObject(b, (void **) &par_b);

   HYPRE_IJVectorAssemble(x);
   HYPRE_IJVectorGetObject(x, (void **) &par_x);


  /*  Print out the system  - files names will be IJ.out.A.XXXXX
       and IJ.out.b.XXXXX, where XXXXX = processor id */
   if (print_system)
   {
      HYPRE_IJMatrixPrint(A, "IJ.out.A");
      HYPRE_IJVectorPrint(b, "IJ.out.b");
   }


   /* Choose a solver and solve the system */

   /* AMG */
   if (globlinsolve == 4)
   {
      int num_iterations;
      double final_res_norm;

      /* Create solver */
      HYPRE_BoomerAMGCreate(&solver);

      /* Set some parameters (See Reference Manual for more parameters) */
      HYPRE_BoomerAMGSetPrintLevel(solver, printLevel);  /* print solve info + parameters */
      HYPRE_BoomerAMGSetNumSweeps(solver, 1);   /* Sweeeps on each level */
      HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
      HYPRE_BoomerAMGSetTol(solver, tollinsolve);      /* conv. tolerance */
      HYPRE_BoomerAMGSetMaxIter(solver,100);
      
      HYPRE_BoomerAMGSetOldDefault(solver); /* Falgout coarsening with modified classical interpolaiton */
      HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* 3 = G-S/Jacobi hybrid relaxation, 0 = jacobi, 1 = Gauss-Seidel, sequential (very slow!) */
      HYPRE_BoomerAMGSetRelaxOrder(solver, 1);   /* uses C/F relaxation */
      HYPRE_BoomerAMGSetCycleType(solver, 1);   /* cycle type. 1=V, 2=W */
      HYPRE_BoomerAMGSetCoarsenType(solver,6);
      HYPRE_BoomerAMGSetStrongThreshold(solver,0.9);

      /* Now setup and solve! */
      HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

      /* Run info - needed logging turned on */
      HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
      HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
      if (myid == 0)
      {
         printf("\n");
         printf("Iterations = %d\n", num_iterations);
         printf("Final Relative Residual Norm = %e\n", final_res_norm);
         printf("\n");
      }

      /* Destroy solver */
      HYPRE_BoomerAMGDestroy(solver);
   }else if(globlinsolve == 5){
      // BICGSTAB with Euclid
      int num_iterations;
      double final_res_norm;

      /* Create solver */
      
      HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD,&solver);
      /* Set some parameters (See Reference Manual for more parameters) */
      HYPRE_ParCSRBiCGSTABSetMaxIter(solver, 10000); /* max iterations */
      HYPRE_ParCSRBiCGSTABSetTol(solver, tollinsolve); /* conv. tolerance */
      HYPRE_ParCSRBiCGSTABSetPrintLevel(solver, printLevel); /* prints out the iteration info */
      HYPRE_ParCSRBiCGSTABSetLogging(solver, 1); /* needed to get run info later */

      // Preconditioner
      HYPRE_EuclidCreate(MPI_COMM_WORLD, &precond);

      HYPRE_ParCSRBiCGSTABSetPrecond(solver, (HYPRE_PtrToParSolverFcn) HYPRE_EuclidSolve, 
         (HYPRE_PtrToParSolverFcn) HYPRE_EuclidSetup, precond);

      /* Now setup and solve! */
      HYPRE_ParCSRBiCGSTABSetup(solver, parcsr_A, par_b, par_x);
      HYPRE_ParCSRBiCGSTABSolve(solver, parcsr_A, par_b, par_x);

      /* Run info - needed logging turned on */
      HYPRE_ParCSRBiCGSTABGetNumIterations(solver, &num_iterations);
      HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm(solver, &final_res_norm);
      if (myid == 0)
      {
         printf("\n");
         printf("Iterations = %d\n", num_iterations);
         printf("Final Relative Residual Norm = %e\n", final_res_norm);
         printf("\n");
      }

      /* Destroy solver */
      HYPRE_ParCSRBiCGSTABDestroy(solver);
      HYPRE_EuclidDestroy(precond);
   }
   else
   {
      if (myid ==0) printf("Invalid solver id specified.\n");
   }
   // copy results
   Copy3DArrayFromHypre(sx, x, grid);

   /* Clean up */
   HYPRE_IJMatrixDestroy(A);
   HYPRE_IJVectorDestroy(b);
   HYPRE_IJVectorDestroy(x);

   /* Finalize MPI*/
   MPI_Finalize();

   return;
}


