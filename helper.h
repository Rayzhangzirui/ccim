#ifndef HELPER_H
#define HELPER_H

#include "grid.h"
#include "matrix.h"
#include <omp.h>
#include <array>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include "global.h"

using std::array;
using std::vector;
using Index = std::array<int,3>;
using Vector3d = std::array<double,3>;

inline Index ToIndex(int* x){
	return Index {x[0],x[1],x[2]};
}


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


// void init_surf_perturb(double ***S, double radius, GridData &grid, int opt, PBData &pb, double tol);

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


vector<vector<double> > randmx(int m, int n);
vector<double> randmx(int n);

double norm1(const vector<vector<double>> &A);

vector<vector<double>> Pointer2dToVector(double **A, int m, int n);
double** VectorToPointer2d(const vector<vector<double>> &A);

//solve Ax = b , dim n, idx 0 to n-1
vector<double> gesolve(const vector<vector<double>> &A, vector<double> b);


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
inline double mean(vector<double> x){
    double sum = 0;
    for (auto& xi : x)
        sum += xi;
    return sum/x.size();
}


inline double variance(vector<double> x){
    double sum = 0;
    double m = mean(x);
    for (auto& xi : x)
        sum += pow(xi - m,2);
    
    return sum/(x.size()-1);
}

inline double rmse(vector<double> x, double y){
    double sum = 0;
    for (auto& xi : x)
        sum += pow(xi - y,2);
    
    return sqrt(sum/x.size());
}

inline double rmse(vector<double> x){
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


inline vector<vector<double>> transpose(const vector<vector<double>> &A){
   int n = A.size();
   vector<vector<double>> B (A);
   for(int i = 0; i < n; i++){
      for(int j = 0; j < n; j++){
         B[i][j] = A[j][i];
      }
   }
   return B;
}

//componentwise
inline vector<double> abs(const vector<double> &x){
   vector<double> v(x);
   for(auto &vi:v){
      vi = abs(vi);
   }
   return v;
}

inline double sum(const vector<double> &x){
   return std::accumulate(x.begin(), x.end(), 0.0);
}


inline vector<double> sign(const vector<double> &x){
   vector<double> v(x);
   for(auto &vi:v){
      vi = (vi>=0)? 1.0: -1.0;
   }
   return v;
}

inline int max_abs_idx(const vector<double> &x){
   vector<double> temp = abs(x);
   return  max_element(temp.begin(),temp.end()) - temp.begin();
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