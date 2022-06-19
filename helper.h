#ifndef HELPER_H
#define HELPER_H

#include "tryvn.h"
#include <omp.h>
#include <array>
#include <vector>
#include <random>

using std::array;
using Index = std::array<int,3>;
using Vector3d = std::array<double,3>;

inline Index ToIndex(int* x){
	return Index {x[0],x[1],x[2]};
}

void cmdLine(int argc, char *argv[]);

//initialize grid
void init_grid(GridData &grid, int nx, double a);

//initialize pb data
void init_pb(PBData &pb, GridData &grid);

// Initialize sparse matrix solver
void init_TempStruct(TempStruct &tmp, GridData &grid);

// initialize PBdata
void init_PBData(PBData &pb, GridData &grid, double epsp = 80.0, double epsm = 1.0);



// subscript to coordinate
inline array<double,3> sub2coord(array<int,3> index, GridData& grid){
	return array<double,3> {grid.a[0]+index[0]*grid.dx[0],
							grid.a[1]+index[1]*grid.dx[1],
							grid.a[2]+index[2]*grid.dx[2] };
}

//distance between two points
inline double dist(array<double,3>x, array<double,3>y){
	return sqrt(pow(x[0]-y[0],2) + pow(x[1]-y[1],2) + pow(x[2]-y[2],2));
}

inline double dist(double *x, double *y){
	return sqrt(pow(x[0]-y[0],2) + pow(x[1]-y[1],2) + pow(x[2]-y[2],2));
}

template<typename T>
inline T evalarray(T*** S, Index ind){
	return S[ind[0]][ind[1]][ind[2]];
}

template<typename T>
inline void setvalarray(T*** S, Index ind, T value){
	S[ind[0]][ind[1]][ind[2]] = value;
}

template<typename T>
inline void addarray(T*** S, Index ind, T value){
	S[ind[0]][ind[1]][ind[2]] += value;
}

inline double GetSign(double*** S, int *x){
	return (evalarray(S,x)>0)? 1.0: -1.0;
}

inline double GetSign(double*** S, Index x){
	return (evalarray(S,x)>0)? 1.0: -1.0;
}


// create sphere surface
void create_sphere(double*** surf, GridData &grid, double radius, array<double,3> center);

//maximum error
double field_max_diff(double ***f1, double ***f2, GridData &grid);


void write_field(std::string fname,double ***matrix, GridData &grid);

//calculate exact u
void get_exact_u(double ***u, double *** surf, GridData &grid);

// print row by col matrix
template<typename T>
void printMat(int row, int col, const T &m){
	for (int i=0; i<row; i++){
		for (int j=0; j<col; j++){
			printf("%.6f ", m[i][j]);
		}
		printf("\n");
	}
}

template<typename T>
void printMat(int row, const T &m){
	for (int i=0; i<row; i++){
		printf("%.6f ", m[i]);
	}
	printf("\n");
}


void init_surf_protein_paper(double ***S, GridData &grid);
void init_surf_protein_dist(double ***S, GridData &grid, int step = 0);

// initialize surface
void init_surf(double ***S, double radius, GridData &grid, int surfopt);

void init_surf_perturb(double ***S, double radius, GridData &grid, int opt, PBData &pb, double tol);

// index at boundary
inline bool atbound(int* ind, GridData &grid){
	for (int d = 0; d < grid.dim; d++){
		if(ind[d]==0||ind[d]==grid.nx[d]){
			return true;
		}
	}
	return false;
}

inline bool atbound(Index ind, GridData &grid){
	return atbound(ind.data(), grid);
}

// index out of obound
inline bool outofbound(int* ind, GridData &grid){
	for (int d = 0; d < grid.dim; d++){
		if(ind[d]<0||ind[d]>grid.nx[d]){
			return true;
		}
	}
	return false;
}

inline bool outofbound(Index ind, GridData &grid){
	return outofbound(ind.data(),grid);
}

inline bool nearinterface(Index ind, double ***S, GridData& grid){
	for(int i : {0,1,2}){
		for(int s : {-1,1}){
			Index nbr = ind;
			nbr[i] = ind[i] + s;
			if(!outofbound(nbr, grid)){ 
				//if nbr not out of bound
				if (S[nbr[0]][nbr[1]][nbr[2]] * S[ind[0]][ind[1]][ind[2]] <= 0.0){
					// if change sign
					return true;
				}
			}
		}
	}
	return false;
}


// print nbr of m around index, m = 1, 3x3x3 nbr, m = 2, 5x5x5 nbr, outside +, inside -
void print_surf(double ***S, int* index, int m);
void print_surf(double ***S, array<int,3> ind, int m);
void print_surf(char ***S, int* index, int m);
void print_surf(char ***S, array<int,3> ind, int m);

void checkwithexactvn(double ***vn, double ***S, PBData &pb, GridData& grid);

void checkcim345DuAll(double ***u, double ***S, PBData &pb, GridData &grid);


int getstatus5debug(double ***S, int *index, GridData &grid);


std::vector<std::vector<double> > randmx(int m, int n);
std::vector<double> randmx(int n);

double norm1(const std::vector<std::vector<double>> &A);
double condest(const std::vector<std::vector<double>> &A);

//solve Ax = b , dim n, idx 0 to n-1
std::vector<double> gesolve(const std::vector<std::vector<double>> &A, std::vector<double> b);

void cim345cond(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, int *index, double ***a, int gamma[][2], double ***S, PBData &pb, GridData &grid);

char yescim5D2All(std::vector<double***> &D2ucoefvec, std::vector<double*> &D2uxcoefvec, std::vector<double*> &D2uxxcoefvec, std::vector<std::vector<int>> &offsetvec, int m, int n, int *index, int mid, double ***S, GridData &grid);

inline void copyMat(double*** toA, double*** fromB, int m, int n, int l){
   for (int i = 0; i <= m; ++i){
      for (int j = 0; j <= n; ++j){
         for (int k = 0; k <= l; ++k){
            toA[i][j][k] = fromB[i][j][k];
         }
      }
   }
}

// for radius statistics
inline double mean(std::vector<double> x){
    double sum = 0;
    for (auto& xi : x)
        sum += xi;
    return sum/x.size();
}


inline double variance(std::vector<double> x){
    double sum = 0;
    double m = mean(x);
    for (auto& xi : x)
        sum += pow(xi - m,2);
    
    return sum/(x.size()-1);
}

inline double rmse(std::vector<double> x, double y){
    double sum = 0;
    for (auto& xi : x)
        sum += pow(xi - y,2);
    
    return sqrt(sum/x.size());
}

inline double rmse(std::vector<double> x){
    double sum = 0;
    for (auto& xi : x)
        sum += pow(xi,2);
    
    return sqrt(sum/x.size());
}


// print plane with + or minus
inline void printPlane(double ***S, GridData &grid, int idx[], int m,int n, int nbr){
	if(m>n){
		std::swap(m,n);
	}
	// row is increasing in n dimension
	// column is increasing in m dimension
	printf("(%i, %i, %i), plane (%i,%i) \n", idx[0], idx[1], idx[2], m, n);
	for(int i = -nbr; i <= nbr; i++){
		for(int j = -nbr; j <= nbr; j++){
			array<int,3> sidx{idx[0],idx[1],idx[2]};
			sidx[m] += i;
			sidx[n] += j;
			if (outofbound(sidx,grid)){
				continue;
			}
			double v = S[sidx[0]][sidx[1]][sidx[2]];
			char m = (v<0.0)?'-':'+';
			std::cout<<m;
		}
		std::cout<<std::endl;
	}
}

inline bool fileExist (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

inline void readMatrix(std::string filename, double ***S, GridData &grid){

	if(!fileExist(filename)){
		std::cerr<<"file do not exit"<<std::endl;
		exit(1);
		return;
	}

	std::ifstream infile(filename);
	//read data from surffile
	for (int i = 0;i <= grid.nx[0];i ++){
		for (int j = 0;j <= grid.nx[1];j ++){
			for (int k = 0;k <= grid.nx[2];k ++){
					infile>>S[i][j][k];
			}
		}
	}

	char next = infile.get();
	if (next==' ' && infile.peek()==EOF){
		std::cout<<"end of file"<<std::endl;
	}else{
		std::cerr<<"not end of file, check grid num."<<std::endl;
		exit(1);
	}
	infile.close();
}

inline bool SameSign(double a, double b){
	return a * b > 0;
}

inline bool SameSide(double*** S, Index a, Index b){
	return SameSign(evalarray(S,a),evalarray(S,b));
}


void OutputVector(std::string filename, double *x, int n);
#endif