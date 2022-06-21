//max difference

#include <string>
#include <sstream>

#include <vector>
#include <algorithm>
#include <numeric>
#include <random>

#include "helper.h"
#include "global.h"
#include "interface.h"
#include "numerics.h"

using namespace std;




//maximum error
double field_max_diff(double ***f1, double ***f2, GridData &grid){
	double max = 0.0;
	double error = max;
	for (int i=0;i<=grid.nx[0];i++){
		for (int j=0;j<=grid.nx[1];j++){
			for (int k=0;k<=grid.nx[2];k++){
				error = abs(f1[i][j][k]-f2[i][j][k]);
				max = (error>max)?error:max;
			}
		}
	}
	return max;
}


void get_exact_u(double ***u, double *** surf, GridData &grid){
   array<double,3> origin = {0,0,0};

   for(int i = 0; i <= grid.nx[0]; i++){
      for(int j = 0; j <= grid.nx[1]; j++){
         for(int k = 0; k <= grid.nx[2]; k++){
            array<double,3> x = sub2coord(array<int,3> {i,j,k}, grid);
            double r = dist(x, origin);
            if(surf[i][j][k] <= 0.0){
               u[i][j][k] = 1.0/(1.0 + pow(r,2));
            }else{
               u[i][j][k] = -1.0/(1.0 + pow(r,2));
            }

         }
      }
   }
}

void write_field(string fname, double ***matrix, GridData &grid){
   ofstream outputfile(fname);
   for (int i=0;i<=grid.nx[0];i++){
     for (int j=0;j<=grid.nx[1];j++){
       for (int k=0;k<=grid.nx[2];k++){
            outputfile << matrix[i][j][k] << " ";
       }
     }
   }
   outputfile.close();
}

// print nbr of m around index, m = 1, 3x3x3 nbr, m = 2, 5x5x5 nbr, outside 1, inside 0
// in right hand coordinant system, view at direction (1,0,0)
void print_surf(double ***S, int* index, int m){
   for (int i = -m; i <= m; i ++){
      printf("i = %d\n", i);
      for (int k = m; k >= -m; k --){
        for (int j = m; j >= -m; j --){
            printf("%c ", (S[index[0]+i][index[1]+j][index[2]+k] < 0)?'-':'+'  );
         }
         printf("\n");
      }
   }
}

void print_surf(double ***S, array<int,3> ind, int m){
  print_surf(S, ind.data(), m);
}

// print nbr of m around index, m = 1, 3x3x3 nbr, m = 2, 5x5x5 nbr, outside 1, inside 0
void print_surf(char ***S, int* index, int m){
   for (int i = -m; i <= m; i ++){
      printf("i = %d\n", i);
      for (int k = m; k >= -m; k --){
        for (int j = m; j >= -m; j --){
            printf("%c ", (S[index[0]+i][index[1]+j][index[2]+k])==0?'-':'+');
         }
         printf("\n");
      }
   }
}

void print_surf(char ***S, array<int,3> ind, int m){
  print_surf(S, ind.data(), m);
}




// void init_surf_perturb(double ***S, double radius, GridData &grid, int opt, PBData &pb, double tol){
//   init_surf(S,  radius,  grid,  opt);
//   if (globperturb == 0)
//     perturb(S,grid.tol,pb,grid);
//   else if (globperturb > 0)
//     perturbstatus(S,grid.tol,globperturb,pb,grid);  
// }


// return BFS, at most n element
vector<Index> BFS(Index idx, double ***S, GridData& grid, int n){
  vector<Index> result;
  vector<Index> queue;
  queue.push_back(idx);
  result.push_back(idx);


  while(!queue.empty() && result.size() < n){
    Index temp = queue.back();
    queue.pop_back();
    for(int r: {0,1,2}){
      for(int s: {-1,1}){
        if(result.size() >= n){
          break;
        }
        Index nbr = temp;
        nbr[r] += s;
        if (SameSide(S, temp, nbr) && find(result.begin(), result.end(), nbr) == result.end()){
          // if same side, not already in result
          queue.push_back(nbr);
          result.push_back(nbr);
        }
      }
    }
  }
  return result;
}
// fillin small void in protein
void FillVoid(double ***S, GridData& grid){
  const int search_max = 7;

  for(int i = 0; i <= grid.nx[0]; i++){
    for(int j = 0; j <= grid.nx[1]; j++){
       for(int k = 0; k <= grid.nx[2]; k++){
        Index idx = {i, j, k};
        
        if(GetSign(S, idx) > 0.0 && nearinterface(idx, S, grid)){
          // if near interface and negative
          vector<Index> connect = BFS(idx, S, grid, search_max);
          if(connect.size() < search_max){
            cout<<"filling void at"<<endl;
            print_surf(S, idx, 3);
            for(Index temp : connect){
              setvalarray(S, temp, -1.0);
            }
          }
        }

      }
    }
  }
}



// void init_surf_protein_dist(double ***S, GridData &grid, int n){
//    vector< vector<double> > p;
//    read_protein(S,grid,p);
//    #pragma omp parallel for collapse(3)
//    for(int i = 0; i <= grid.nx[0]; i++){
//       for(int j = 0; j <= grid.nx[1]; j++){
//          for(int k = 0; k <= grid.nx[2]; k++){
//             array<double,3> x = sub2coord(array<int,3> {i,j,k}, grid);
//             S[i][j][k] = sqrt(3.0) * 2 * fabs(grid.a[0]);

//             for(int n = 0; n < p.size(); n++){
//                array<double,3> px = {p[n][0],p[n][1],p[n][2]};//px is coordinate
//                double d =  dist(x,px) - p[n][3];
//                S[i][j][k] = min(S[i][j][k],d);
//             }
            
//          }
//       }
//    }

//    TempStruct tmp;
//    init_TempStruct(tmp, grid);

//    for(int i = 0; i<n; i++){
//       advanceheat(S, tmp, grid);
//    }
// }


//create random matrix
vector<vector<double>> randmx(int m, int n){
   vector<vector<double>> A(m, vector<double>(n,0.0));
   static default_random_engine e;
   static uniform_real_distribution<double> u(-1,1);
   for(int i = 0; i<m; i++){
      for(int j = 0; j<n; j++){
         A[i][j] = u(e);
      }
   }
   return A;
}

vector<double> randmx(int n){
   vector<double> x(n, 0.0);
   static default_random_engine e;
   static uniform_real_distribution<double> u(-1,1);
   for(int i = 0; i<n; i++){
         x[i] = u(e);
   }
   return x;
}

//calculate 1 norm, of matrix A dim n+1 x n+1, max column sum
double norm1(const vector<vector<double>> &A){
   double max = 0.0;
   int m = A.size();
   int n = A[0].size();
   for(int j = 0; j < n; j++){
      double temp = 0.0;
      for(int i = 0; i < m; i++){
         temp += abs(A[i][j]);
      }
      max = (temp>max)?temp:max;
   }
   return max;
}

// 
vector<vector<double>> Pointer2dToVector(double **A, int m, int n){
   vector<vector<double>> B(m, vector<double>(n,0.0));
   for(int i = 0; i < m; i++){
      for(int j = 0; j < n; j++){
         B[i][j] = A[i][j];
      }
   }
   return B;
}


double** VectorToPointer2d(const vector<vector<double>> &A){
   int m = A.size();
   int n = A[0].size();
   double **B = matrix(m-1,n-1);
   for(int i = 0; i < m; i++){
      for(int j = 0; j < n; j++){
         B[i][j] = A[i][j];
      }
   }
   return B;
}

//solve Ax = b , dim n, idx 0 to n-1
vector<double> gesolve(const vector<vector<double>> &A, vector<double> b){
   int n = b.size();
   int PLR[n], PLC[n];
   double **LU = matrix(n-1,n-1);
   double **B =  VectorToPointer2d(A);
   
   vector<double> x(n,0.0);
   
   gecp0(LU,PLR,PLC,B,n-1,n-1);
   
   forwardbacksub0(x.data(),b.data(),LU,PLR,PLC,n-1);
   
   free_matrix(LU,n-1,n-1);
   free_matrix(B,n-1,n-1);
   return x;
}





void OutputVector(string filename, double *x, int n)
{
   int i, r;
   ofstream outfile(filename);
   outfile.precision(16);
   for (i = 0; i < n; i++){
      outfile << i  << " " << scientific << x[i] << endl; 
   }
         
}

