#ifndef ICIM_H
#define ICIM_H

#include "helper.h"
#include <algorithm>
#include <array>

using std::min;
using std::max;
using std::abs;
using std::array;
using std::vector;
using std::cout;
using std::endl;


// unit increment in dim 
inline array<int,3> UnitIncrement(Index ind, int dim, int inc){
	Index e = ind;
	e[dim] += inc;
	return e;
}


template<typename T>
void print(T* g){
	printf("(%i,%i,%i)", g[0],g[1],g[2]);
}

template<typename T>
void print(array<T,3> g){
	printf("(%i,%i,%i)", g[0],g[1],g[2]);
}


template<typename T>
inline T Dot(array<T,3> a, array<T,3> b){
	return (a[0]*b[0])+(a[1]*b[1])+(a[2]*b[2]);
}

template<typename T>
inline array<T,3> Add(array<T,3> a, array<T,3> b){
	return array<T,3> {a[0]+b[0] , a[1]+b[1] , a[2]+b[2] };
}

template<typename T>
inline array<T,3> Minus(array<T,3> a, array<T,3> b){
	return array<T,3> {a[0]-b[0] , a[1]-b[1] , a[2]-b[2] };
}

template<typename T>
inline double Norm2(array<T,3> a){
	return sqrt(Dot(a,a));
}

inline double NormInf(array<double,3> a){
	return max(max(abs(a[0]),abs(a[1])),abs(a[2]));
}

inline array<double,3> Mult(double a, array<double,3> b){
	return array<double,3> {a*b[0] , a*b[1] , a*b[2] };
}

inline array<double,3> Div(array<double,3> b, double a){
	return array<double,3> {b[0]/a , b[1]/a , b[2]/a };
}



inline bool IsGhost(int ind[], double*** sign_surf, double*** S){
	return !SameSign(evalarray(sign_surf,ind),evalarray(S,ind));
}

inline bool IsGhost(Index ind, double*** sign_surf, double*** S){
	return !SameSign(evalarray(sign_surf,ind),evalarray(S,ind));
}

inline bool isnan(Vector3d v){
	return std::any_of(v.cbegin(), v.cend(), [](double e) { return isnan(e); });
}




template <typename T>
inline void SetMatrix(T ***S, T value, int c, int r, int f){
	if(S){
		for(int i = 0 ; i <= c; i++){
			for(int j = 0 ; j <= r; j++){
				for(int k = 0 ; k <= f; k++){
					S[i][j][k] = value;
				}
			}
		}	
	}
}

template <typename T>
inline void SetField(T ***S, GridData &grid, T value){
	SetMatrix(S, value, grid.nx[0], grid.nx[1], grid.nx[2]);
}

template <typename T>
inline void SetMatrix(T **S, T value, int c, int r){
	if(S){
		for(int i = 0 ; i <= c; i++){
			for(int j = 0 ; j <= r; j++){
					S[i][j] = value;
			}
		}	
	}
}

template <typename T>
inline void SetMatrix(T *S, T value, int c){
	if(S){
		for(int i = 0 ; i <= c; i++){
			S[i] = value;
		}	
	}
}


// initialize sign surface cS from  level set function dS
inline void SignSurf(double*** signS, double *** dS, GridData &grid){
	for(int i = 0 ; i <= grid.nx[0]; i++){
		for(int j = 0 ; j <= grid.nx[1]; j++){
			for(int k = 0 ; k <= grid.nx[2]; k++){
				signS[i][j][k] = (dS[i][j][k] < 0.0)? -1.0: 1.0;
			}
		}
	}
}

inline void PrintCoef(double*** u, int N){
	for( int iloc = 0; iloc <= N; iloc ++){
		for( int jloc = 0; jloc <= N; jloc ++){
			for( int kloc = 0; kloc <= N; kloc ++){
				if (abs(u[iloc][jloc][kloc])>1e-9){
					printf("%.6f  v(%i,%i,%i)\n", u[iloc][jloc][kloc],iloc,jloc,kloc);
				}
			}
		}
	}
}

inline void PrintCoef(double* u, int N){
	for( int iloc = 0; iloc <= N; iloc ++){
		if (abs(u[iloc])>1e-9){
			printf("%.6f  v(%i)\n", u[iloc], iloc);
		}
	}
}

// get direct nbr excluding it self
inline std::vector< array<int,3> > GetDirectNeighbors( Index ind, GridData& grid){
	std::vector< array<int,3> > v;
	v.reserve(6);
	for (int r = 0; r < 3; r++){
		for (int s = -1; s <= 1; s +=2){
			Index t = UnitIncrement(ind, r, s);
			if(!outofbound(t,grid)){
				v.push_back(t);
			}
		}
	}
	return v;
}

// get cube nbr
inline std::vector< array<int,3> > GetCubeNeighbors( Index ind, GridData& grid, int side){
	std::vector< array<int,3> > v;
	
	for(int i = -side; i <=side; i ++){
		for(int j = -side; j <=side; j ++){
			for(int k = -side; k <=side; k ++){
				Index t = Add(ind,Index{i,j,k});
				if(!outofbound(t,grid)){
					v.push_back(t);
				}
			}
		}
	}
	return v;
}


// print nbr of m around index, m = 1, 3x3x3 nbr, m = 2, 5x5x5 nbr, outside 1, inside 0
inline void PrintFlipSurf(double *** flip_surf, double *** original_surf, int* index, int m){
   for (int i = -m; i <= m; i ++){
      printf("i = %d\n", i);
      for (int k = m; k >= -m; k --){
        for (int j = m; j >= -m; j --){
        	Index temp = Add(ToIndex(index), Index{i,j,k});
        	if(IsGhost(temp, flip_surf, original_surf)){
        		printf(" (%s) ", (flip_surf[temp[0]][temp[1]][temp[2]])<0.0? "-":"+");	
        	}else{
        		printf("  %s  ", (flip_surf[temp[0]][temp[1]][temp[2]])<0.0? "-":"+");	
        	}
            
         }
         printf("\n");
      }
   }
}

array<int,3> Indicator(Index ind, double*** signS, GridData &grid);

// check if is type 1 exception 
bool istype1(Index ind, double*** signS, GridData &grid);

// check if is interior point
bool isinterior(Index ind, double*** signS, GridData &grid);

// check if is type 2 exception
bool istype2(Index ind, double*** signS, GridData &grid);

// type 0 for interior, type 1 and 2 as defined in paper, eqn(12), 3 for regular cim, 4 for boundary
int GetType( Index ind, double*** signS, GridData &grid );

void flip(double*** sign_surf, double*** original_surf, PBData &pb, GridData &grid, int maxIter, int depth, bool cond = false);

void linearsystem_icim(SparseElt2**** &A, double ***b, double ***S, double*** sign_surf, PBData &pb, GridData &grid);

bool IsUsualCim2(Index ind, double***S, GridData& grid);

bool HasUsualMixD2u(double*** u, Index ind, int k, int l, double*** S, GridData& grid);

bool HasCoupleMixD2u(double*** u, double* uxxcoef, Index ind, int k, int l, double*** S, GridData& grid);

bool HasShiftMixD2u(double*** u, Index ind, int k, int l, double*** S, GridData& grid);

bool HasMixD2u(double ***u,  double *uxxcoef, int m, int n, Index index, double ***S, GridData &grid);

// use shift to approximate Du
bool HasShiftDuAndMix(double ***ucoef, double* uxxcoef, Index ind, int dur, double*** S, GridData &grid);

vector<Index> GetShiftD2uNbr(Index ind, int dur, double*** S, GridData &grid);

bool IsIcim2(Index ind, double***S, GridData&grid);

void GetIcimDu(double ****ucoef, double **uxxcoef, double *** S, int *index, int rstar, int sstar, double alpha, GridData &grid);

void CheckIcimDu(double*** u_ghost, double *** u, double*** sign_surf, double ***S,  PBData &pb, GridData &grid);

void GetIcimDuGridpoint(double ****ucoef, double **uxxcoef, double *** S, int *index, GridData &grid);

void ExpandCoeff(double coef, Index a, double *** A, Index b, double*** B, int N);

// reconstruct uxx at grid point
array<double,3> ComputeD2uLookup(Index index, double*** S, GridData& grid);

// reconstruct mix D2u at grid point
void ComputeMixD2u(double** D2u, double*** u_ghost, double *** sign_surf, double *** S, Index index, GridData &grid);

// reconstruct Du at grid point
void ComputeDuGridpoint(double Du[], double*** u_ghost,  double *** sign_surf, double *** S, Index index, GridData &grid);

// reconstruct u at interface
double ComputeIcimU(double*** u_ghost, double *** sign_surf, double *** S, Index index, int rstar, int sstar, double alpha, GridData &grid);

// construct coupling equation, calculte determinant
double DetCouplingMatrix(Index index, double *** flipping_S, double *** S, PBData&pb, GridData &grid);

// reconstruct u at grid point
void reconstruct(double*** psi_true, double*** S, double*** psi_ghost, double*** sign_surf, GridData& grid);

// check u err at grid point
void CheckErrGrid(double *** u, double ***S,  PBData &pb, GridData &grid);

// check u err at interface
void CheckIcimU(double*** u_ghost, double *** u, double*** sign_surf, double ***S,  PBData &pb, GridData &grid);

void GetCouplingMatrix(double**M, double*d, double ****f, int *index, double ***S, double*** sign_surf, PBData &pb, GridData &grid);

void ComputeByExtrapolate(double& uval, Vector3d& du, double*** u_ghost, double *** sign_surf, double *** S, Index index, int rstar, int sstar, double alpha, GridData &grid);

void CheckErrGeneral(double*** psi_true, double*** S, double*** psi_ghost, double*** sign_surf, PBData &pb, GridData &grid);

template<typename T>
inline void PrintAarry(T *x, int n){
	cout<<"(";
	for(int i = 0; i < n - 1; i++){
		cout<<std::fixed<<std::setprecision(6)<<x[i]<<", ";
	}
	cout<<x[n-1]<<")";
}

// print error for n element vector
inline void PrintError(int n, const char* var1, double* d1, const char* var2, double* d2, const char* var3){
	vector<double> err(n, 0.0);
	for(int i = 0; i < n; i++){
		err[i] = abs(d1[i]-d2[i]);
	}

	cout<<var1<<" = ";
	PrintAarry(d1, n);
	cout<<" ";

	cout<<var2<<" = ";
	PrintAarry(d2, n);
	cout<<" ";

	cout<<var3<<" = ";
	PrintAarry(err.data(), n);
	cout<<" ";

	cout<<endl;
}

// print error for scalar
inline void PrintError(const char* var1, double d1, const char* var2, double d2, const char* var3){
	cout<<std::fixed<<std::setprecision(6)<<var1<<" = "<<d1<<" ";
	cout<<std::fixed<<std::setprecision(6)<<var2<<" = "<<d2<<" ";
	cout<<std::fixed<<std::setprecision(6)<<var3<<" = "<<abs(d1-d2)<<" ";
	cout<<endl;
}
#endif