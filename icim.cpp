#include "icim.h"
#include "tryvn.h"
#include "interface.h"
#include "global.h"
#include <unordered_map>
#include <list>
#include <iterator>
using namespace std;

#define DEBUG

int CUBENBR = 2;

struct Info{
	double d[3];
	double f[3][5][5][5];
	double M[3][3];
	double d2u[3];
	double sign;
};


Info SaveInfo(double** M, double* d, double**** f){
	Info info;
	for(int r = 0; r < 3; r++){
		info.d[r] = d[r];
		for(int i = 0 ; i <= 4; i++){
			for(int j = 0 ; j <= 4; j++){
				for(int k = 0 ; k <= 4; k++){
					info.f[r][i][j][k] = f[r][i][j][k];
				}
			}
		}
	}

	for(int i = 0 ; i <= 2; i++){
			for(int j = 0 ; j <= 2; j++){
				info.M[i][j] = M[i][j];
		}
	}

	return info;
};

// Retreve info
void GetInfo(double** M, double* d, double**** f, const Info& info){
	
	for(int r = 0; r < 3; r++){
		d[r] = info.d[r];
		for(int i = 0 ; i <= 4; i++){
			for(int j = 0 ; j <= 4; j++){
				for(int k = 0 ; k <= 4; k++){
					f[r][i][j][k] = info.f[r][i][j][k];
				}
			}
		}
	}

	for(int i = 0 ; i <= 2; i++){
		for(int j = 0 ; j <= 2; j++){
			M[i][j] = info.M[i][j];
		}
	}
};

unordered_map<int, Info> D2uMap;


// equation (12) in ICIM paper, return ind[r] = 3 if the r-dim touch boundary
array<int,3> Indicator(Index ind, double*** sign_surf, GridData &grid){
	
	array<int,3> g = {0,0,0};
	for(int r = 0; r < 3; r ++){

		if (outofbound(UnitIncrement(ind,r,1), grid) || outofbound(UnitIncrement(ind,r,-1), grid) ){
			// if outofbound, 3
				g[r] = 3;
		}else{
			// if not out of bound
			if( SameSign(evalarray(sign_surf,UnitIncrement(ind,r,1)),evalarray(sign_surf,ind)) 
				&& SameSign(evalarray(sign_surf,UnitIncrement(ind,r,-1)),evalarray(sign_surf,ind)) ){
				g[r] = 2;
			}else if( !SameSign(evalarray(sign_surf,UnitIncrement(ind,r,1)),evalarray(sign_surf,ind)) &&
			          !SameSign(evalarray(sign_surf,UnitIncrement(ind,r,-1)),evalarray(sign_surf,ind)) ){
				g[r] = 0;
			}else{
				g[r] = 1;
			}

		}		
	}
	return g;
}


// check if is type 1 exception 
bool istype1(Index ind, double *** sign_surf, GridData &grid){

	for( int k = 0; k < 3; k++){
		Index tind = ind; // temporary index
		array<int,3> g0 = Indicator(ind, sign_surf, grid);

		tind[k] = ind[k]-1;
		array<int,3> gm = Indicator(tind, sign_surf, grid); //g(x_{i-e_k})

		tind[k] = ind[k]+1;
		array<int,3> gp = Indicator(tind, sign_surf, grid); //g(x_{i+e_k})

		// type 1 exception
		if (g0[k]==0 || gm[k]==0 || gp[k] == 0){
			return true;
		}
	}
	return false;
}

// check if is interior point
bool isinterior(Index ind, double*** S, GridData &grid){
	array<int,3> g = Indicator(ind, S, grid);
	return (g[0]==2 && g[1]==2 && g[2]==2);
}

// check if is type 2 exception
bool istype2(Index ind, double*** S, GridData &grid){
	if (istype1(ind, S, grid)){
		return false;
	}

	for (int m = 0; m < grid.dim; m++){
		for (int n = m+1; n < grid.dim; n++){
			int sk2[4];

			char has_sk2 = yessk2(sk2, m, n, ind.data(), S, grid);
			
			// 
			// if(ind[0]==eindex[0] && ind[1]==eindex[1] && ind[2]==eindex[2]){
			// 	printf("in %i %i plane, sk2=%i\n",m,n,has_sk2);
			// 	printPlane(S, grid, ind.data(), m, n, 2);
			// }
			// 

			if (!has_sk2){
				return true;
			}
		}
	}
	return false;
}

// type 0 for interior, type 1 and 2 as defined in paper, eqn(12), 3 for regular cim, 4 for boundary
int GetType( Index ind, double ***S, GridData &grid ){

	if(atbound(ind, grid)){
		return 4;
	}

	if (isinterior(ind, S, grid)){
		return 0;
	}

	if(istype1(ind, S, grid)){
		return 1;
	}

	if(istype2(ind, S, grid)){
		return 2;
	}

	return 3;
}


void AtomFlip(Index ind, double*** S, GridData& grid){
	S[ind[0]][ind[1]][ind[2]] = -1 * S[ind[0]][ind[1]][ind[2]];

	vector<Index> nbrs = GetDirectNeighbors(ind, grid);

	// reset do-not-flip flag for nbrs
	for(auto nbr : nbrs){
		if( abs(S[nbr[0]][nbr[1]][nbr[2]])>1.5){
			S[nbr[0]][nbr[1]][nbr[2]] = S[nbr[0]][nbr[1]][nbr[2]]/2.0;
		}
	}

	// debug info
	// 
	// int type = GetType(ind, S, grid);
	// printf("flip at (%i,%i,%i)  type %i\n", ind[0],ind[1],ind[2], type);
	// print_surf(S,ind,2);
	// 
}

// flip ind and its negbhors upto certain recursion depth
// if if flipbydet is true, use determinant as threshold for flipping
// if marked as do-not-flip, set value to +/- 2
void LocalFlip(Index ind, double *** S, double *** flipping_S,  PBData&pb, GridData &grid, int depth, int &counter, bool flipbydet){
	if (depth == 0){
		return;
	}

	vector<Index> nbrs = GetDirectNeighbors(ind, grid);
	int type = GetType( ind, flipping_S, grid); // use flip S
	if (type == 4){
		// skip boundary

	}else if (type == 1){
		//if type 1, see table 1
		
		array<int,3> g = Indicator(ind, flipping_S, grid);
		Index gsort = g;
		std::sort(gsort.begin(),gsort.end());
		
		
		if (gsort == array<int,3>{0,0,0}){ 
			//flip
			AtomFlip(ind, flipping_S, grid);
			++counter;
			LocalFlip(ind,  S, flipping_S, pb, grid, --depth, counter, flipbydet);
		
		}if(gsort == array<int,3>{0,0,1}){ // estimate determinant and flip

			if (abs(flipping_S[ind[0]][ind[1]][ind[2]])>1.5){
				// if marked as do-not-flip, condinue
				return;
			}
			
			bool fliphere = true;
			
			AtomFlip(ind, flipping_S, grid); // first time flip
			
			if(flipbydet){
				// calculate condition number
				double det = DetCouplingMatrix(ind, flipping_S, S, pb, grid);

				if(det < 0.1265){
					printf("det = %f at (%i,%i,%i) Do not flip.\n",det, ind[0],ind[1],ind[2] );
					AtomFlip(ind, flipping_S, grid); // flip back
					fliphere = false;
					flipping_S[ind[0]][ind[1]][ind[2]] = 2 * flipping_S[ind[0]][ind[1]][ind[2]]; // mark as -2 or 2

				}
			}



			if(fliphere){
				// if it is actually fliped, flip nbr
				++counter;
				LocalFlip(ind,  S, flipping_S, pb, grid, --depth, counter, flipbydet);		
			}
			
			
			


		}else if(gsort == array<int,3>{0,1,1} ){

			Index dummy_idx;
			int dummy_r;
			int r = std::distance(g.begin(), std::find(g.begin(), g.end(), 0)); //find with dimension is 0
			
			// -- first shift then flip
			vector<Index> nbr_for_d2u = GetShiftD2uNbr(ind, r, flipping_S, grid);
			if( nbr_for_d2u.empty()){
				AtomFlip(ind, flipping_S, grid);
				++counter;
				LocalFlip(ind,  S, flipping_S, pb, grid, --depth, counter, flipbydet);
			}
			
			// -- flip no matter what
			// bool fliphere = true;
			// AtomFlip(ind, flipping_S, grid); // first time flip
			// if(flipbydet){
			// 	// calculate condition number
			// 	double det = DetCouplingMatrix(ind, flipping_S, S, pb, grid);

			// 	if(det < 0.1265){
			// 		printf("det = %f at (%i,%i,%i) Do not flip.\n",det, ind[0],ind[1],ind[2] );
			// 		AtomFlip(ind, flipping_S, grid); // flip back
			// 		fliphere = false;
			// 	}
			// }
			// if(fliphere){
			// 	// if it is actually fliped, flip nbr
			// 	++counter;
			// 	LocalFlip(ind,  S, flipping_S, pb, grid, --depth, counter, flipbydet);		
			// }

		}


		//for dim r where g[r] == 0, if can not find shift, then have to flip
		for(int r = 0; r < 3; r ++){
			if(g[r] == 0){
				vector<Index> nbr_shift_d2u = GetShiftD2uNbr(ind, r, flipping_S, grid);
				if( nbr_shift_d2u.empty()){
					printf("g = (%i,%i,%i) at (%i,%i,%i).  No shift for D2u[%i] , has to flip.\n", g[0], g[1], g[2], ind[0], ind[1], ind[2], r);
					AtomFlip(ind, flipping_S, grid);
					++counter;
					LocalFlip(ind,  S, flipping_S, pb, grid, --depth, counter, flipbydet);
				}			
			}
		}

	}else if(type == 2){
		// if type 2, see figure 6
		// if case 4 5 6, either more than 2 nbr(up down left right) has G = 0, flip

		for(int k = 0; k < 3; k++){
			for(int l = k+1; l < 3; l++){
				// iterater over plane-k,l
				int count_plane = 0;

				array<int,2> dim = {k,l}; 
				array<int,2> dir = {-1,1};
				vector< Index > index_to_flip; // list of index to flip
				for(int r: dim){
					for(int s: dir){
						Index nbr_idx = UnitIncrement(ind, r, s); // neighbor in k l plane
						array<int,3> G = Indicator(nbr_idx, S, grid);

						int other_dim = (r == dim[0])? dim[1] : dim[0];
						
						if (G[other_dim] == 0){
							count_plane++;
							index_to_flip.push_back( nbr_idx );
						}
					}
				}

				if (count_plane>=2){
					
					for( auto idx: index_to_flip){
						AtomFlip(idx, flipping_S, grid);
						++counter;
						LocalFlip(idx, S, flipping_S, pb, grid, --depth, counter, flipbydet);
					}

					AtomFlip(ind, flipping_S, grid); // flip current
					++counter;
					LocalFlip(ind, S,  flipping_S, pb, grid, --depth, counter, flipbydet);
				}				 
			}
		} // end of kl loop
	}// end of type 2



	return;
}


// method to flip iteratively
void flip(double*** sign_surf, double*** original_surf, PBData &pb, GridData &grid, int maxIter, int depth, bool cond){
	if (cond){
		cout<<"Flipping with coupling matrix determinant"<<endl;
	}else{
		cout<<"Simple flipping without considering determinant"<<endl;
	}
	int total_flip = 0;
	int iter_num = maxIter; // number of iteration

	for(int iter = 0; iter < maxIter; iter++){
		int count_iter = 0; // counter for flip action

		for(int i = 0 ; i <= grid.nx[0]; i++){
			for(int j = 0 ; j <= grid.nx[1]; j++){
				for(int k = 0 ; k <= grid.nx[2]; k++){
					Index ind = {i,j,k};
					
					LocalFlip(ind, original_surf, sign_surf,  pb, grid, depth, count_iter, cond);	
					
					
				}
			}
		}//end of ijk loop

		total_flip += count_iter;

		if (count_iter == 0){
			// if no flip in this iteration, break out of max iter
			iter_num = iter;
			break;
			
		}
	}


	int total_sign_change = 0;
	for(int i = 0 ; i <= grid.nx[0]; i++){
			for(int j = 0 ; j <= grid.nx[1]; j++){
				for(int k = 0 ; k <= grid.nx[2]; k++){
					if (!SameSign(sign_surf[i][j][k], original_surf[i][j][k])){
						total_sign_change ++;
					}
				}
			}
		}//end of ijk loop

	cout<< "Flip " << total_flip << " in "<< iter_num << " iterations"<< endl;
	cout<< "Total sign change = " << total_sign_change <<endl;
	return;
}


bool HasTwoNeighborBothSide(Index ind, int dim, int s, double*** S){
	double a[4];
	a[0] = evalarray(S, UnitIncrement(ind,dim,-s));//1 step back
	a[1] = evalarray(S, ind);// here
	a[2] = evalarray(S, UnitIncrement(ind,dim,s));//1 step forward
	a[3] = evalarray(S, UnitIncrement(ind,dim,2*s));// 2step forward
	return SameSign(a[0],a[1]) && !SameSign(a[1],a[2]) && SameSign(a[2],a[3]);
}

bool HasTwoNeighborWithGap(Index ind, int dim, int s, double*** S){
	double a[4];
	a[0] = evalarray(S, UnitIncrement(ind,dim,-s));//1 step back
	a[1] = evalarray(S, ind);// here
	a[2] = evalarray(S, UnitIncrement(ind,dim,s));//1 step forward
	a[3] = evalarray(S, UnitIncrement(ind,dim,2*s));// 2step forward
	return SameSign(a[0],a[1]) && !SameSign(a[1],a[2]) && SameSign(a[1],a[3]);
}

bool HasTwoNeighbor(Index ind, int dim, double*** S){
	double a[3];
	int s = 1;
	a[0] = evalarray(S, UnitIncrement(ind,dim,-s));//1 step back
	a[1] = evalarray(S, ind);// here
	a[2] = evalarray(S, UnitIncrement(ind,dim,s));//1 step forward
	return SameSign(a[0],a[1]) && SameSign(a[1],a[2]);
}


// approx second derivative by it nbr, used for coupling equation
// return a vector of Index that can be used to approximate D2u
vector<Index> GetShiftD2uNbr(Index ind, int dur, double*** S, GridData &grid){
	vector<int> dim = {0,1,2};
	vector<Index> result;
	vector<Index> nbr_order;

	vector<pair<double,Index> > dist_ind_pairs;
	
	Index range{CUBENBR,CUBENBR,CUBENBR}; // only look at cube nbr
	range[dur] = CUBENBR - 1; // along dur range can only be 1, if 2, then may step out of 5-point stencil
	for( int i = -range[0]; i <= range[0]; i ++){
		for( int j = -range[1]; j <= range[1]; j ++){
			for( int k = -range[2]; k <= range[2]; k ++){
				Index nbr = Add(ind, Index {i,j,k});
				if (SameSide(S, ind, nbr) && HasTwoNeighbor(nbr, dur, S)){
					// same side with two nbr
					double dist = Norm2(Minus(nbr, ind));
					dist_ind_pairs.push_back(make_pair(dist, nbr));
				}
			}
		}
	}

	// sort by distance
	sort(dist_ind_pairs.begin(),dist_ind_pairs.end());

	for(auto e : dist_ind_pairs){
		result.push_back(e.second);
	}



	return result;
}


// struct FD{
// 	vector<double> coef;
// 	vector<Index> stencil;
// 	Index center;
// };

// std::vector<FD>  GetUsualMixD2u(double*** u,  Index ind, int k, int l, double*** S, GridData& grid){
// 	vector<FD> fds;
// 	int combimation[9][4] = {
// 		{-1, 1, -1, 1},
// 		{-1, 1, -1, 0},
// 		{-1, 1,  0, 1},
// 		{-1, 0, -1, 1},
// 		{ 0, 1, -1, 1},
// 		{-1, 0, -1, 0},
// 		{-1, 0,  0, 1},
// 		{ 0, 1, -1, 0},
// 		{ 0, 1,  0, 1} };
// 	for( int i = 0; i < 9; i++){

// 	}


// }


// has usual mixed derivative, four point squre stencile 
bool HasUsualMixD2u(double*** u,  Index ind, int k, int l, double*** S, GridData& grid){
	int mid = 2;
	int t, N = 2*mid, sindex[grid.dim], rindex[grid.dim], nindex[grid.dim], sk2[4];
	
	
	
	if (yessk2(sk2,k,l,ind.data(),S,grid))//if D2u in (m,n) plane at index has approximation in same side
	{
		for (int t = 0; t < grid.dim; t++)
		{
			sindex[t] = mid;
			nindex[t] = N;
		}

		if(u){
			// if pointer is not null, write
			SetMatrix(u, 0.0, N, N, N);
			getD2(u,k,l,sk2,sindex,nindex,grid);// approximate D2u in terms of u-value only
		}

		return true;
	}
	return false;
}


// get O(h) apprx of Du from nbr points
bool HasShiftMixD2u(double*** ucoef, Index ind, int k, int l, double*** S,  GridData& grid){
	int mid = 2;
	int N = 2 * mid;
	SetMatrix(ucoef, 0.0, N, N, N);
	
	bool found = false;

	array<int,3> dim = {0,1,2};	
	vector<Index> nbr_order;
	auto it = find_if(dim.begin(), dim.end(), [k,l](int v)->bool{return v!=k &&v!=l;});
	int o = *it; // out of pane dimension

	Index range = {mid - 1, mid - 1, mid - 1};
	range[o] = mid;
	
	// push back same side, order by distance, prefer out of plane
	for(int i = -range[0]; i <= range[0]; i ++ ){
		for(int j = -range[1]; j <= range[1]; j ++ ){
			for(int k = -range[2]; k <= range[2]; k ++ ){
				Index nbr_idx = Add(ind, Index{i,j,k});
				if(SameSide(S, ind, nbr_idx)){
					nbr_order.push_back(nbr_idx);
				}
			}
		}
	}
	auto comp = [ind,o](Index &a,  Index &b){
		double da = Norm2(Minus(a,ind));
		double db = Norm2(Minus(b,ind));
		if(abs(da-db) > 1e-9){
			// return smaller 2 norm
			return da < db ;	
		}else{
			// tie breaker, return out of plane
			return abs(a[o] - b[o]) > 0;
		}
	};
	std::sort( nbr_order.begin(), nbr_order.end(), comp);  



	double*** temp = matrix(N,N,N);
	for(Index nbr_idx : nbr_order ){
		if( !found && HasUsualMixD2u(temp, nbr_idx, k, l, S, grid)){
			// shift result by grid in dim r, dir s
			ExpandCoeff(1.0, ind, ucoef, nbr_idx, temp,  N);
			found = true;
		}

		if(found){
			break;
		}
	}	
	
	free_matrix(temp, N, N, N);
	return found;
}




bool HasCoupleMixD2u(double*** u, double* uxxcoef, Index ind, int k, int l, double*** S, GridData& grid){
	int mid = 2;
	int N = 2 * mid;

	SetMatrix(u, 0.0, N, N, N);
	SetMatrix(uxxcoef, 0.0, 2);

	double hsqured = (grid.dx[k]*grid.dx[l]);

	for ( int sk = -1; sk <=1; sk +=2){
		for ( int sl = -1; sl <=1; sl +=2){
			//restart every time

			SetMatrix(u, 0.0, N, N, N);
			SetMatrix(uxxcoef, 0.0, 2);

			Index center = {mid,mid,mid};
			
			if ( SameSide(S, ind, UnitIncrement(ind,k,sk)) && SameSide(S, ind, UnitIncrement(ind,l,sl)) ){
				Index local_stepk = UnitIncrement(center,k,sk); // 1 step in k
				Index local_stepl = UnitIncrement(center,l,sl); // 1 step in l
				Index local_stepkl = UnitIncrement(local_stepk,l,-sl); // 1 step in k, then opp step in l
				Index local_steplk = UnitIncrement(local_stepl,k,-sk); // 1 step in l, then opp step in k 
				
				if (SameSide(S, ind, Add(ind, Minus(local_steplk,center)))) {
					// couple with ukk
					uxxcoef[k] = (sk*sl) ;

					setvalarray(u, center, 1.0/hsqured/(sk*sl) );
					setvalarray(u, local_stepk,  -1.0/hsqured/(sk*sl) );
					setvalarray(u, local_stepl,  1.0/hsqured/(sk*sl) );
					setvalarray(u, local_steplk,  -1.0/hsqured/(sk*sl) );
					return true;
				}

				if (SameSide(S, ind, Add(ind, Minus(local_stepkl,center)))) {
					uxxcoef[l] = (sk*sl);
					
					setvalarray(u, center, 1.0/hsqured/(sk*sl) );
					setvalarray(u, local_stepk,  1.0/hsqured/(sk*sl) );
					setvalarray(u, local_stepl,  -1.0/hsqured/(sk*sl) );
					setvalarray(u, local_stepkl,  -1.0/hsqured/(sk*sl) );
					return true;
				}				
			}
		}
	}
	return false;
}

bool IsUsualCim2(Index ind, double***S, GridData&grid){
	array<int,3> g = Indicator(ind, S, grid);

	// for each dimension
	for(int r = 0; r < 3; r++){

		if(g[r] == 3){
			cerr<<"pass boundary points to IsUsualCim2"<<endl;
			exit(1);
		}

		// has two nbr
		if(g[r] == 2){
			continue;
		}

		// has interface both side
		if(g[r] == 0){
			return false;
		}

		// has interface on one side
		if(g[r] == 1){
			int s=1;
			if(SameSign(evalarray(S,ind), evalarray(S,UnitIncrement(ind, r, s)))){
				s = -1;
			} 

			if (!HasTwoNeighborBothSide(ind, r, s, S)){
				return false;
			}
			
		}		


		// find mixed derivative
		for(int m = 0; m < 3; m++){
			for(int n = m + 1; n < 3; n++){
				if (!HasUsualMixD2u( nullptr, ind, m, n, S, grid)){
					return false;
				}
			}
		}


	}

	return true;
}



bool IsIcim2(Index ind, double***S, GridData&grid){
	array<int,3> g = Indicator(ind, S, grid);
	
	for(int r = 0; r < 3; r++){
		if(g[r] == 3){
			cerr<<"pass boundary points to IsUsualCim2"<<endl;
			exit(1);
		}
	}


	if (IsUsualCim2(ind, S, grid)){
		return true;
	}else{
		int type = GetType(ind, S, grid);


		if(type==1){
			bool has_fix = true;

			for(int r = 0; r < 3; r++){
				// each dim should has a fix

				if(g[r] == 0){
					vector<Index> nbr_shift_d2u = GetShiftD2uNbr( ind, r, S, grid);
					if (nbr_shift_d2u.empty()){
						// if both side has interface, but not shift, then fail
						return false;
					}
				}else if (g[r]==1){
					// if one side has interface, 
					// interface side
					int s = (SameSide(S, ind, UnitIncrement(ind, r, -1) ))? 1 : -1;

					if ( ! (HasTwoNeighborBothSide(ind, r, s, S) || HasTwoNeighborWithGap(ind, r, s, S))){
						return false;
					}
				}else{
					// both side no interface, not a problem
				}
			}

		}else if (type==2){
			// each plane should have a fix;
			for(int m = 0; m < 3; m ++){
				for(int n = m + 1; n < 3; n ++){
					
					int dummy_s, dummy_r;

					if (!HasMixD2u(nullptr, nullptr, m, n, ind,  S, grid)) {
						// if do not have couple mix deri or usual mix Du
						return false;
					}
				}
			}
		}
	}

	return true;
}

// get mixed derivitive
bool HasMixD2u(double ***u,  double *uxxcoef, int m, int n, Index index, double ***S, GridData &grid)
{
	if (m>n){
		swap(m, n);
	}
	int mid = 2;
	int t, N = 2*mid, sindex[grid.dim], rindex[grid.dim], nindex[grid.dim], sk2[4];
	double thesign;
	if(u){
		SetMatrix(u, 0.0, N, N, N);	
	}
	
	if(uxxcoef){
		SetMatrix(uxxcoef, 0.0, 2);	
	}
	
	

	if (HasUsualMixD2u(u, index, m, n, S, grid)){
		return true;
	
	}else if (HasCoupleMixD2u(u, uxxcoef, index, m, n, S, grid)){
		return true;

	}else if(HasShiftMixD2u(u, index, m, n, S, grid)){
		return true;

	}else{
		fprintf(stderr, "[HasMixD2u] no D2u[%i][%i] at (%i,%i,%i)\n", m, n, index[0],index[1],index[2]);
		print_surf(S, index, 3);		
	}
	return false;

}

// A is centered at a, B is centered at b, transfer coeff of B to A, assume both size N x N x N
// coeff of b should be copy to b-a
void ExpandCoeff(double coef, Index a, double *** A, Index b, double*** B, int N){

	Index offset = Minus(b,a);
	for( int i = 0; i <= N; i ++){
		for( int j = 0; j <= N; j ++){
			for( int k = 0; k <= N; k ++){
				Index Bind = {i,j,k}; // Bind is an index in B
				Index Aind = Add(Bind , offset); //Aind is Bind's local coord in A 
				if (all_of(Aind.begin(), Aind.end(), [N](int e) { return e<=N && e>=0; })){
					// if Aind is in bound
					A[Aind[0]][Aind[1]][Aind[2]] += coef * B[Bind[0]][Bind[1]][Bind[2]];
				}else{
					if (abs(evalarray(B,Bind))>1e-9){
						cerr<<"loosing infomation when Expanding Coeff"<<endl;
						exit(1);
					}
				}
			}
		}
	}
}


// when approximating Du[dur]
// xnbr = x + sr h e_r
// D[dur](x) = D[dur](xnbr)  - sr h D[dur][r](x) + O(h^2)
// D[dur](xnbr) = D^0 (xnbr) + O(h^2) 
// 			 or = D^s (xnbr) + 0.5 s h D[dur][dur](x) + O(h^2) , s = 1 backward diff, s = -1 forward diff
bool HasShiftDuAndMix(double ***ucoef, double* uxxcoef, Index ind, int dur, double*** S, GridData &grid){
	
	int mid = 2;
	int N = mid * 2;
	Index center = {mid, mid, mid};
	SetMatrix(ucoef, 0.0, N, N, N);
	SetMatrix(uxxcoef, 0.0, 3);


	double*** mix_d2u_ucoef = matrix(N, N, N);
	double mix_d2u_uxxcoef[3];


	bool found = false;
	// find nbr with central diff
	for(int r: {0, 1, 2}){
		if(r == dur){
			continue;
		}

		if( !found && HasMixD2u(mix_d2u_ucoef, mix_d2u_uxxcoef, r, dur, ind, S, grid)){

			for(int s : {-1, 1}){
				Index nbr = UnitIncrement(ind, r, s);
				Index nbr_loc = UnitIncrement(center, r, s);
				if(SameSide(S, ind, nbr)){
					Index g = Indicator(nbr, S, grid);
					if(g[dur] == 2){
						// can do central diff
						addarray(ucoef, UnitIncrement(nbr_loc, dur , 1), 1.0 / (2.0 * grid.dx[dur]));
						addarray(ucoef, UnitIncrement(nbr_loc, dur , -1), -1.0 / (2.0 * grid.dx[dur]));

						found = true;
					}else if(g[dur] == 1){
						// can do forward backward diff
						int snbr = SameSide(S, UnitIncrement(nbr,dur,1), nbr)? -1: 1;// which side is interface

						addarray(ucoef, nbr_loc, snbr / grid.dx[dur]);
						addarray(ucoef, UnitIncrement(nbr_loc, dur , -snbr), - snbr / grid.dx[dur]);
						uxxcoef[dur] += 0.5 * snbr * grid.dx[dur];

						found = true;
					}else{
						// 2 interface, go to next nbr
						continue;
					}
				}

				if(found){
					ExpandCoeff(-s * grid.dx[dur], ind, ucoef, ind, mix_d2u_ucoef, N);
					break;
				}

			}	

		}
		
	}


	free_matrix(mix_d2u_ucoef, N, N, N);

	if(ind[0]==eindex[0] && ind[1]==eindex[1] && ind[2]==eindex[2]){
		printf("[HasShiftDuAndMix] (%i,%i,%i) \n",ind[0], ind[1], ind[2]);
		// for each component of Du
		double dummy_D1uxcoef[3] = {0};
		double thesign = GetSign(S,ind);
		double approx = evalcoef(0.0, ucoef,dummy_D1uxcoef,uxxcoef,ind.data(),0,0,0.0,mid,S,grid);
		double exact_Du = getDu(ind.data(),dur,0,0,0,thesign,grid);
		double err = approx - exact_Du;
		printf("grid Du[%i] = %f, exact = %f, loc trunc err =  %f\n", dur, approx, exact_Du, err);
			
	}

	return found;
}


// Approximate Du in dur dimension
// a nbr in rr-dim has central diff in dur and mix Du(rr,r), used to find Du at grid point
// xnbr = x + sr h e_r+ sr2 h e_r2 + sr3 h e_r3
// D[r](x) = D[r](xnbr)  - sr h D[r][r](x) - sr2 h D[r][r2](x) - sr3 h D[r][r3](x) + O(h^2)
// D[r](xnbr) = D^0 (xnbr) + O(h^2) 
// or D^s (xnbr) + 0.5 s h D[r][r](x) + O(h^2) , s = 1 backward diff, s = -1 forward diff

bool HasShiftDuAndMixGeneral(double ***ucoef, double* uxxcoef, Index ind, int dur, double*** S, GridData &grid){
	
	int mid = 2;
	int N = mid * 2;
	Index center = {mid, mid, mid};
	SetMatrix(ucoef, 0.0, N, N, N);
	SetMatrix(uxxcoef, 0.0, 3);


	double*** mix_d2u_ucoef = matrix(N, N, N);
	double mix_d2u_uxxcoef[3];


	bool found = false;


	vector<Index> nbr_order;
	Index range = {mid , mid , mid};
	range[dur] = mid - 1; // smaller range in du dimension

	// push back same side, order by distance, prefer central difference
	for(int i = -range[0]; i <= range[0]; i ++ ){
		for(int j = -range[1]; j <= range[1]; j ++ ){
			for(int k = -range[2]; k <= range[2]; k ++ ){
				Index nbr_idx = Add(ind, Index{i,j,k});
				Index g = Indicator(nbr_idx, S, grid);

				if(SameSide(S, ind, nbr_idx) && g[dur]!=0){
					// same side and has neighbor
					nbr_order.push_back(nbr_idx);
				}
			}
		}
	}
	auto comp = [&](Index &a,  Index &b){
		double da = Norm2(Minus(a,ind));
		double db = Norm2(Minus(b,ind));
		if(abs(da-db) > 1e-9){
			// return smaller 2 norm
			return da < db ;	
		}else{
			// prefer central diff
			return HasTwoNeighbor(a, dur, S);
		}
	};
	std::sort( nbr_order.begin(), nbr_order.end(), comp);  

	if(nbr_order.empty()){
		fprintf(stderr, "no nbr candidates at (%i,%i,%i)\n", ind[0], ind[1], ind[2]);
		print_surf(S, ind, mid);
		exit(1);
	}

	// approximate second derivative
	
	Index nbr_to_eval;
	for(Index nbr: nbr_order){
		
		SetMatrix(ucoef, 0.0, N, N, N);
		SetMatrix(uxxcoef, 0.0, 3);

		int count_valid = 0; // count if each dimension sucess;
		Index s = Minus(nbr, ind);
		for(int r: {0, 1, 2}){
			if(s[r] == 0){
				count_valid ++;
				continue;
			}else{
				if(r == dur){
					mix_d2u_uxxcoef[r] += -s[r];
					count_valid ++;
				}else{
					if(HasMixD2u(mix_d2u_ucoef, mix_d2u_uxxcoef, r, dur, ind, S, grid)){
						count_valid ++;
						ExpandCoeff(-s[r] * grid.dx[dur], ind, ucoef, ind, mix_d2u_ucoef, N);
					}
				}

			}
		}
		if (count_valid == 3){
			nbr_to_eval = nbr;
			found = true;
			break;
		}
	}

	if(!found){
		fprintf(stderr, "[HasShiftDuAndMixGeneral] no D2u at (%i,%i,%i)\n", ind[0], ind[1], ind[2]);
		print_surf(S, ind, mid);
		exit(1);
	}

	// approximated Du[dur] at nbr
	Index nbr_loc = Add(center, Minus(nbr_to_eval, ind));
	if(HasTwoNeighbor(nbr_to_eval, dur, S)){
		addarray(ucoef, UnitIncrement(nbr_loc, dur , 1), 1.0 / (2.0 * grid.dx[dur]));
		addarray(ucoef, UnitIncrement(nbr_loc, dur , -1), -1.0 / (2.0 * grid.dx[dur]));

	}else{

		int snbr = SameSide(S, UnitIncrement(nbr_to_eval,dur,1), nbr_to_eval)? -1: 1;// which side is interface
		addarray(ucoef, nbr_loc, snbr / grid.dx[dur]);
		addarray(ucoef, UnitIncrement(nbr_loc, dur , -snbr), - snbr / grid.dx[dur]);
		uxxcoef[dur] += 0.5 * snbr * grid.dx[dur];
	}

	// find nbr with central diff

	free_matrix(mix_d2u_ucoef, N, N, N);

	if(ind[0]==eindex[0] && ind[1]==eindex[1] && ind[2]==eindex[2]){
		printf("[HasShiftDuAndMixGeneral] (%i,%i,%i) \n",ind[0], ind[1], ind[2]);
		// for each component of Du
		double dummy_D1uxcoef[3] = {0};
		double thesign = GetSign(S,ind);
		double approx = evalcoef(0.0, ucoef,dummy_D1uxcoef,uxxcoef,ind.data(),0,0,0.0,mid,S,grid);
		double exact_Du = getDu(ind.data(),dur,0,0,0,thesign,grid);
		double err = approx - exact_Du;
		printf("grid Du[%i] = %f, exact = %f, loc trunc err =  %f\n", dur, approx, exact_Du, err);
			
	}

	return found;
}

// Approx Du at grid piont Du[r] = ucoef[r][:][:][:] u[:] + uxxcoef[r][:] uxx[:]
void GetIcimDuGridpoint(double ****ucoef, double **uxxcoef, double *** S, int *index, GridData &grid){

	int mid = 2;
	int N = mid * 2;
	Index center = {mid, mid, mid};
	
	for(int r = 0; r < 3; r++){
		SetMatrix(ucoef[r], 0.0, N, N, N);	
		SetMatrix(uxxcoef[r], 0.0, 2);
	}
	
	

	Index ind = ToIndex(index);

	Index g = Indicator(ind, S, grid);
	

	for(int r = 0; r < 3; r ++){
		// for each Du[r]

		if (g[r] == 2){
			// has two nbr along r
			// approx Du_r(x) by central diff
			setvalarray(ucoef[r], UnitIncrement(center, r,  1), 1.0 / (2.0 * grid.dx[r]));
			setvalarray(ucoef[r], UnitIncrement(center, r, -1), - 1.0 / (2.0 * grid.dx[r]));
		}else if (g[r] == 1){
			// has interface in r direction
			// approx Du_r(x) by forward/backward difference, 
			// Du[r](xhat) = (1/h) D[r,s] u(x) + s (h/2) D2_rr (x) + alpha h D2_{r,r*} (x)
			int s = SameSide(S, ind, UnitIncrement(ind, r, 1))? -1 : 1; 

			setvalarray(ucoef[r], center, s/grid.dx[r]);
			setvalarray(ucoef[r], UnitIncrement(center, r,  -s),  -s/grid.dx[r]);

			uxxcoef[r][r] = 0.5 * s * grid.dx[r];

		}else if(g[r] == 0){
			
			//D[r](x)  + sr h D[r][r](x) + srr h D[r][rr](x) = D[r](x + sr h e_r+ srr h e_sr) 
			if(HasShiftDuAndMix(ucoef[r], uxxcoef[r], ind, r, S, grid)){

			}else if(HasShiftDuAndMixGeneral(ucoef[r], uxxcoef[r], ind, r, S, grid)){

			}else{
				fprintf(stderr,"No shift and mixDu available to approx Du[%i] at grid point (%i, %i, %i)\n", r, ind[0], ind[1], ind[2]);
				print_surf(S, ind, 2);
				exit(1);
			}
		}
	}

	
	if(ind[0]==eindex[0] && ind[1]==eindex[1] && ind[2]==eindex[2]){
		printf("[GetIcimDuGridpoint] (%i,%i,%i) \n",index[0], index[1], index[2]);
		// for each component of Du
		double dummy_D1uxcoef[3] = {0};
		for(int m = 0; m < 3; m++){
			double thesign = GetSign(S,ind);
			double approx = evalcoef(0.0, ucoef[m],dummy_D1uxcoef,uxxcoef[m],ind.data(),0,0,0.0,mid,S,grid);
			cout<<"ucoef"<<endl;
			PrintCoef(ucoef[m],N);
			cout<<"uxxcoef"<<endl;
			PrintCoef(uxxcoef[m],2);
			double exact_Du = getDu(ind.data(),m,0,0,0,thesign,grid);
			double err = approx - exact_Du;
			printf("Grid Du: g = %i, Du[%i] = %f, exact = %f, loc trunc err =  %f\n", g[m], m, approx, exact_Du, err);
		}	
	}	
}



// Get O(h^2) approxmination of Du at interface 
void GetIcimDu(double ****ucoef, double **uxxcoef, double *** S, int *index, int rstar, int sstar, double alpha, GridData &grid){
	int mid = 2;
	int N = mid * 2;

	Index ind = ToIndex(index);

	Index g = Indicator(ind, S, grid);

	// approx cim Du at grid point
	
	GetIcimDuGridpoint(ucoef, uxxcoef, S, index, grid);
	
	for(int r = 0; r < 3; r ++){
		if(r == rstar){
			uxxcoef[r][r] += alpha * sstar * grid.dx[rstar];
		}else{
			double coef_mixDu =  alpha * sstar *  grid.dx[rstar];

			double *** mixDu_ucoef = matrix(N, N, N); // mixDu in r, rstar plane in terms of u
			double mixDu_uxxcoef[3] = {0}; // mixDu in r, rstar plane in terms of uxx
			if (HasMixD2u(mixDu_ucoef,  mixDu_uxxcoef, r, rstar, ind , S, grid)){
				
				ExpandCoeff(coef_mixDu, ind, ucoef[r], ind, mixDu_ucoef, N);

				for (int i = 0; i < 3; i ++){
					uxxcoef[r][i] += coef_mixDu * mixDu_uxxcoef[i];
				}
			}else{
				fprintf(stderr, "cannot approximate mixDu at plane (%i,%i) at (%i, %i, %i), rstar = %i, sstar = %i\n",
				 r, rstar, index[0], index[1], index[2], rstar, sstar);
				exit(1);
			}

			free_matrix(mixDu_ucoef, N, N, N);
		}	
	}

	
	if(ind[0]==eindex[0] && ind[1]==eindex[1] && ind[2]==eindex[2]){
		printf("[GetIcimDu] (%i,%i,%i) rstar = %i, star = %i\n",index[0], index[1], index[2], rstar, sstar);
		// for each component of Du
		double dummy_D1uxcoef[3] = {0};
		for(int m = 0; m < 3; m++){
			double thesign = GetSign(S,ind);
			double approx = evalcoef(0.0, ucoef[m],dummy_D1uxcoef,uxxcoef[m],ind.data(),0,0,0.0,mid,S,grid);
			double exact_Du = getDu(ind.data(),m,rstar,sstar,alpha,thesign,grid);
			double err = approx - exact_Du;
			printf("Du[%i] = %f, exact = %f, loc trunc err =  %f\n", m, approx, exact_Du, err);
		}	
	}

}



void GetCouplingMatrix(double**M, double*d, double ****f, int *index, double ***S, double*** sign_surf, PBData &pb, GridData &grid)
{
	int r, s, i, j, k, tmps, t, m, mid = 2, N = 2*mid;
	int rindex[grid.dim], tindex[grid.dim], Narray[grid.dim];
	for (r = 0; r < grid.dim; r++)
		Narray[r] = N;
	
	int PLR[grid.dim], PLC[grid.dim];
	double alpha[grid.dim], beta, tangent[grid.dim], normal[grid.dim];
	double bk, rthere, rhere, ethere, ehere, ehat, C;
	double temp[grid.dim], a[4];
	int sk[3] = {0,0,0};
	double sigma, thesign, tau, Dtau[grid.dim];

	// look at flipped sign to classify inside or outside
	if (evalarray(sign_surf,index) < 0)
	{
		ehere = pb.epsilonm;
		ethere = pb.epsilonp;
		thesign = -1.0;
	}
	else
	{
		ehere = pb.epsilonp;
		ethere = pb.epsilonm;
		thesign = 1.0;
	}
	

	Index g = Indicator(ToIndex(index),sign_surf,grid);
	
	
	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		cout<<endl<<"[GetCouplingMatrix]"<<endl;
		cout<<"original surface"<<endl;
		print_surf(S, index, 2);
		cout<<"flipped surface"<<endl;
		print_surf(sign_surf, index, 2);
	}
	

	for (r = 0; r < grid.dim; r++)
		rindex[r] = index[r];
	for (r = 0; r < grid.dim; r++)
		tindex[r] = mid;
	for (r = 0; r < grid.dim; r++) // dim by dim discretization
	{
		// one interface
		if (g[r] == 1){

			sk[r] = (SameSide(sign_surf, ToIndex(index), UnitIncrement( ToIndex(index), r, -1) ))? 1 : -1; // direction of flip interface
			
			if (HasTwoNeighborBothSide(ToIndex(index), r, sk[r], sign_surf)){
				// can use cim2 in this dimension	
				
				//--icim, if the other side is flipped, different alpha
				rindex[r] = index[r]+sk[r];

				if (IsGhost(rindex,sign_surf,S) && !IsGhost(index,sign_surf,S)){
					// the other side is ghost, this side is not, interface is 1 step forward away
					Index here = ToIndex(index);
					Index there = ToIndex(rindex);
					int step = 0;
					while (IsGhost(there,sign_surf,S)){
						here[r] += sk[r];
						there[r] += sk[r];
						step++;
					}
					getinterfaceinfo(alpha[r],tangent,normal,S,here.data(),there.data(),grid);
					alpha[r] += step;

					if(step>1){
						printf("%i consectutive flipping at (%i,%i,%i) along r = %i, s = %i\n", step, index[0],index[1],index[2],r,sk[r]);
					}
					
				}else if ( IsGhost(index,sign_surf,S) && !IsGhost(rindex,sign_surf,S)){
					// this side is ghost, the other side is not, interface is behind, shift back
					Index here = ToIndex(index);
					Index there = ToIndex(rindex);
					int step = 0;
					while (IsGhost(here,sign_surf,S)){
						here[r] -= sk[r];
						there[r] -= sk[r];
						step++;
					}
					getinterfaceinfo(alpha[r],tangent,normal,S,here.data(),there.data(),grid);
					alpha[r] = -step + alpha[r];

					if(step>1){
						printf("%i consectutive flipping at (%i,%i,%i) along r = %i, s = %i\n", step, index[0],index[1],index[2],r,sk[r]);
					}
				}else if(IsGhost(index,sign_surf,S) && IsGhost(rindex,sign_surf,S)){
					fprintf(stderr, "Both side are ghost at (%i,%i,%i) along r = %i, s = %i\n", index[0],index[1],index[2],r,sk[r]);
					exit(1);
				}
				else{
					// both not ghost
					getinterfaceinfo(alpha[r],tangent,normal,S,index,rindex,grid);
				}
				rindex[r] = index[r];
				//--
				

				beta = 1.0-alpha[r];
				ehat = (beta+beta*beta)*(0.5+alpha[r])*ehere+
						 (alpha[r]+alpha[r]*alpha[r])*(0.5+beta)*ethere;
				rhere = ehere/ehat;
				rthere = ethere/ehat;
				a[0] = (beta+beta*beta)*rhere+alpha[r]*(1.0+2.0*beta)*rthere; //a_{i,-1} to a_{i,2} eqn(10) cim paper
				a[1] = -(beta+beta*beta)*rhere-(1.0+alpha[r])*(1.0+2.0*beta)*rthere;
				a[2] = (1.0+beta)*(1.0+beta)*rthere;
				a[3] = -beta*beta*rthere;
				bk = -(beta+beta*beta)*(rthere-rhere);
				C = bk*tangent[r]; // C = b_k (t_k . e_k), eqn(34), coeff for grau(u-)
				
				//L_k operator
				tindex[0] = mid; tindex[1] = mid; tindex[2] = mid;
				for (int s = 0; s < 4; s++)
				{
					tindex[r] = mid+(s-1)*sk[r];
					setvalarray(f[r],tindex,evalarray(f[r],tindex)+a[s]); //u coeff from L operator eqn(32)
				}
				
				//J_k operator
				getsigma(sigma,index,r,sk[r],alpha[r],normal,pb,grid);
				gettau(tau,index,r,sk[r],alpha[r],grid);
				getDtau(Dtau,index,r,sk[r],alpha[r],grid);
				d[r] = thesign*(abs(sk[r])*(1.0+2.0*beta)*rthere*tau+
						sk[r]*(beta+beta*beta)*grid.dx[r]*
						(sigma/ehat*normal[r]+
						rthere*getdotprod(Dtau,tangent,grid.dim)*tangent[r]));// J_k in eqn(36) cim paper

				
				if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
					printf("[GetCouplingMatrix] (%i,%i,%i) g[%i]==1, s = %i, usual CIM2\n",index[0],index[1],index[2],r,sk[r]);
					printf("Dim %i, dir %i, cim2 , alpha = %f \n", r, sk[r], alpha[r]);
					double L = 0;
					for( int iloc = 0; iloc <= N; iloc ++){
						for( int jloc = 0; jloc <= N; jloc ++){
							for( int kloc = 0; kloc <= N; kloc ++){
								if (abs(f[r][iloc][jloc][kloc])>1e-6){
									int globindex[3] = {index[0] + iloc - mid, index[1] + jloc - mid, index[2] + kloc - mid};
									L += f[r][iloc][jloc][kloc] * getu(globindex, 0, 0, 0, GetSign(sign_surf,globindex), grid);
								}
							}
						}
					}

					printf("alpha= %f, tangent = [%f,%f,%f]\n", alpha[r], tangent[0],tangent[1],tangent[2]);
					
					double Du_dot_t = 0;

					for(int s = 0; s < 3; s ++){
						Du_dot_t += getDu(index, s, r, sk[r], alpha[r], thesign, grid) * tangent[s];

					}

					double apprxD2u = 1.0 / pow(grid.dx[r],2) * (d[r] + L + sk[r] * C * grid.dx[r] * Du_dot_t); // use cim2 formular with exact Du-
					double exactD2u = getD2u(index,r,r,0,0,0.0,thesign,grid);

					printf("exact D2u = %f, apprx D2u = %f, loc truc err = %f （with exact Du)\n", exactD2u, apprxD2u, abs(exactD2u - apprxD2u));
				}
				

				// dot( grad(u-(xhat)), t)

				double**** D1ucoef = new double ***[grid.dim]; 
				for (int r = 0; r < grid.dim; r++)
					D1ucoef[r] = matrix(N,N,N);
				double** D1uxxcoef = matrix(grid.dim-1,grid.dim-1);

				GetIcimDu(D1ucoef ,D1uxxcoef, sign_surf, index, r, sk[r], alpha[r], grid);

				for(int j = 0; j < 3; j ++){
					 
					for( int iloc = 0; iloc <= N; iloc ++){
						for( int jloc = 0; jloc <= N; jloc ++){
							for( int kloc = 0; kloc <= N; kloc ++){
								tindex[0] = iloc;
								tindex[1] = jloc;
								tindex[2] = kloc;
								setvalarray(f[r], tindex, evalarray(f[r],tindex) +  (C * tangent[j] * sk[r] * grid.dx[r]) * D1ucoef[j][iloc][jloc][kloc] );
							}
						}
					}
				}

				M[r][r] = 1;
				for(int j = 0; j < 3; j ++){ // from  dot product 
					for(int i = 0; i < 3; i ++){ // for each D2u_ii
						M[r][i] = M[r][i] -  1/ pow(grid.dx[r],2) * (C * tangent[j] * sk[r] * grid.dx[r]) * D1uxxcoef[j][i];
					}
					
				}


				if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
					
					double tempuxcoef[3] = {0};
					double tempuxxcoef[3] = {0}; 
					double lhs = 0;

					double rhs = evalcoef(d[r],f[r], tempuxcoef, tempuxxcoef, index, 0, 0, 0, mid, sign_surf, grid );
					rhs *= 1.0 / pow(grid.dx[r],2);

					for(int j = 0; j < 3; j++){
						lhs += M[r][j] * getD2u(index,j,j,0,0,0.0,thesign,grid);
					}
					printf("lhs = %f, rhs = %f, residual = %f\n", lhs, rhs, abs(rhs-lhs));
				}
				
				// free temp
				free_matrix(D1uxxcoef, 2, 2);
				for (int r = 0; r < grid.dim; r++)
					free_matrix(D1ucoef[r], N,N,N);

				


					
				



			} else if (HasTwoNeighborWithGap(ToIndex(index), r, sk[r], sign_surf)){
				// s=-1, s=0, s= 2 are same side
				tindex[0] = mid; tindex[1] = mid; tindex[2] = mid;//reset tindex
				sk[r] = (SameSide(sign_surf, ToIndex(index), UnitIncrement( ToIndex(index), r, -1) ))? 1 : -1; // direction of interface
						
				tindex[r] = mid;
				setvalarray(f[r],tindex,evalarray(f[r],tindex) - 1.0); //icim eqn(20), coef of u_i

				tindex[r] = mid - sk[r];
				setvalarray(f[r],tindex,evalarray(f[r],tindex) + 2.0/3.0); // icim eqn(20), coef of u_{i-ei}
				
				tindex[r] = mid + 2 * sk[r] ;
				setvalarray(f[r],tindex,evalarray(f[r],tindex) + 1.0/3.0); // icim eqn(20), coef of u_{i + 2ei}
				
				M[r][r] = 1.0;

				if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
					printf("[icim2all] (%i,%i,%i) g[%i]==1, central differencing with gap s = %i\n",index[0],index[1],index[2],r,sk[r]);

					double tempuxcoef[3] = {0};
					double tempuxxcoef[3] = {0}; 
					double lhs = 0;

					double rhs = evalcoef(d[r],f[r], tempuxcoef, tempuxxcoef, index, 0, 0, 0, mid, sign_surf, grid );
					rhs *= 1.0 / pow(grid.dx[r],2);

					for(int j = 0; j < 3; j++){
						lhs += M[r][j] * getD2u(index,j,j,0,0,0.0,thesign,grid);
					}
					printf("lhs = %f, rhs = %f, residual = %f\n", lhs, rhs, abs(rhs-lhs));
				}
				

			} else{
				cerr<<"no fix for uxx g = 1"<<endl;
				exit(1);
			}
			// end of case g[r] == 1 
		} else if(g[r] == 0){
			// two interface on dim r, shift to nbr point
			Index shift_nbr;
			vector<Index> shift_nbr_avail = GetShiftD2uNbr(ToIndex(index), r ,sign_surf, grid);
			if (!shift_nbr_avail.empty()){
				shift_nbr = shift_nbr_avail[0];
				Index mididx = {mid,mid,mid};
				Index temp = Add(mididx,Minus(shift_nbr,ToIndex(index)));//local index of shift_nbr centered at index

				setvalarray(f[r],temp,evalarray(f[r],temp) - 2.0); //icim eqn(20), coef of u_i - sk 
				
				setvalarray(f[r],UnitIncrement(temp,r,1),evalarray(f[r],UnitIncrement(temp,r,1)) + 1.0); //icim eqn(20), coef of u_i - sk + sj

				setvalarray(f[r],UnitIncrement(temp,r,-1),evalarray(f[r],UnitIncrement(temp,r,-1)) + 1.0); //icim eqn(20), coef of u_i - sk - sj

				M[r][r] = 1.0;


				
				if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
					printf("[GetCouplingMatrix] (%i,%i,%i) g[%i]==0, apprx D2u at (%i,%i,%i)\n",
						index[0],index[1],index[2], r, shift_nbr[0],shift_nbr[1],shift_nbr[2]);

					double tempuxcoef[3] = {0};
					double tempuxxcoef[3] = {0}; 
					double lhs = 0;

					double rhs = evalcoef(d[r],f[r], tempuxcoef, tempuxxcoef, index, 0, 0, 0, mid, sign_surf, grid );
					rhs *= 1.0 / pow(grid.dx[r],2);

					for(int j = 0; j < 3; j++){
						lhs += M[r][j] * getD2u(index,j,j,0,0,0.0,thesign,grid);
					}
					printf("lhs = %f, rhs = %f, residual = %f\n", lhs, rhs, abs(rhs-lhs));
				}
				

			}else{
				fprintf(stderr, "no fix for uxx, g[%i] = 0 at (%i,%i,%i)",r,index[0],index[1],index[2]);
				print_surf(sign_surf, ToIndex(index), 2);
				exit(1);
			}

			// end of case g[r] == 0
		} else {
			// g[r]==2, central differencing
			tindex[0] = mid; tindex[1] = mid; tindex[2] = mid;

			M[r][r] = 1.0;
			for (int s = -1; s <= 1; s += 2)
			{
				tindex[r] = mid+s;
				setvalarray(f[r],tindex,evalarray(f[r],tindex)+1.0);
			}
			tindex[r] = mid;
			setvalarray(f[r],tindex,evalarray(f[r],tindex)-2.0);

		}
	
	}

}

void icim2all(SparseElt2**** &A, double ***b, int *index, double ***S, double*** sign_surf, PBData &pb, GridData &grid)
{
	int r, s, i, j, k, tmps, t, m, mid = 2, N = 2*mid;
	int rindex[grid.dim], tindex[grid.dim], Narray[grid.dim];
	for (r = 0; r < grid.dim; r++)
		Narray[r] = N;
	double ****f;
	double **LU, **M, value, Sval;
	int PLR[grid.dim], PLC[grid.dim];
	double alpha[grid.dim], beta, tangent[grid.dim], normal[grid.dim];
	double bk, rthere, rhere, ethere, ehere, ehat, C;
	double temp[grid.dim], a[4];
	int sk[3] = {0,0,0};
	double sigma, thesign, tau, Dtau[grid.dim];
	

	LU = matrix(grid.dim-1,grid.dim-1);
	M = matrix(grid.dim-1,grid.dim-1);
	f = matrix(grid.dim-1,N,N,N); // rhs of coupling linear system,  r-the row M urr =  d[r] + f[r][idx] u[idx]
	double d[3] = {0,0,0};
	
	if (evalarray(sign_surf,index) < 0)
	{
		ehere = pb.epsilonm;
		ethere = pb.epsilonp;
		thesign = -1.0;
	}
	else
	{
		ehere = pb.epsilonp;
		ethere = pb.epsilonm;
		thesign = 1.0;
	}
	


	GetCouplingMatrix(M, d, f, index, S, sign_surf, pb, grid);
	

	gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);
	
	double det = abs(LU[0][0]*LU[1][1]*LU[2][2]);
	// store D2u information
	Info info = SaveInfo(M, d, f);
	info.sign = thesign;
	D2uMap[ sub2ind(index,grid.nx,grid.dim) ] = info;
	

	
	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		printf("[icim2all] (%i,%i,%i)\n",index[0],index[1],index[2]);
		double v[3] = {0};
		double tempuxcoef[3] = {0};
		double tempuxxcoef[3] = {0};
		for(int j = 0; j < 3; j++){
			printf("d[%i] = %.6f\n ", j, d[j]);
			PrintCoef(f[j],N);
		  //double evalcoef(double u0, double ***ucoef, double *uxcoef, double *uxxcoef, int *index, int rstar, int sstar, double alpha, int mid, double thesign, GridData grid)
		  v[j] = 1/pow(grid.dx[0],2) * evalcoef(d[j],f[j], tempuxcoef, tempuxxcoef, index, 0, 0, 0, mid, sign_surf, grid );
		}

		
		cout<<"M matrix det = "<< det <<endl;
		printMat(3,3,M);

		cout<<"exact rhs"<<endl;
		printMat(3,v);
		forwardbacksub0(v,v,LU,PLR,PLC,grid.dim-1);
		

		double exactD2u[3] = {0};
		for(int j = 0; j < 3; j++){
		  exactD2u[j] = getD2u(index,j,j,0,0,0.0,thesign,grid);  
		}
		cout<<"solved,       exact,       error"<<endl;
		for(int j = 0; j < 3; j++){
		  cout<<std::setw(12)<< v[j]<<", "<<std::setw(12)<<exactD2u[j]<<", "<<std::setw(12)<<v[j] - exactD2u[j]<<endl;
		}
		cout <<" -eps lap(u) - f = " <<- ehere * (v[0]+v[1]+v[2])-getf(index,0,0,0.0,thesign,pb,grid)<<endl;
	}
	// end of local truncation error
  


	// Set up A and b
	for (i = 0; i < grid.dim; i++)
		tindex[i] = 0;
	while (tindex[0] <= N)
	{
		for (r = 0; r < grid.dim; r++)
			temp[r] = evalarray(f[r],tindex);
		forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
		value = 0.0;
		for (r = 0; r < grid.dim; r++)
			value += -ehere*temp[r]/(grid.dx[r]*grid.dx[r]);
		if (value != 0.0)
		{
			for (r = 0; r < grid.dim; r++)
				rindex[r] = index[r]-mid+tindex[r];
			sparse2(index,rindex,A,value,grid);
		}

		for (r = 0; r < grid.dim; r++)
			setvalarray(f[r],tindex,temp[r]);

		(tindex[grid.dim-1])++;
		for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
		{
			tindex[i] = 0;
			(tindex[i-1])++;
		}
	}


	forwardbacksub0(d,d,LU,PLR,PLC,grid.dim-1);
	value = 0.0;
	for (r = 0; r < grid.dim; r++)
		value += ehere*d[r]/(grid.dx[r]*grid.dx[r]);
	setvalarray(b,index,evalarray(b,index)+value);




	free_matrix(M,grid.dim-1,grid.dim-1);
	free_matrix(LU,grid.dim-1,grid.dim-1);
	free_matrix(f,grid.dim-1,4,4,4);
	
}

double DetCouplingMatrix(Index index, double *** flipping_S, double *** S, PBData&pb, GridData &grid){
	int N = 4;
	int PLR[3],PLC[3];
	double** LU = matrix(grid.dim-1,grid.dim-1);
	double**M = matrix(grid.dim-1,grid.dim-1);
	double**** f = matrix(grid.dim-1,N,N,N); // rhs of coupling linear system,  r-the row M urr =  d[r] + f[r][idx] u[idx]
	double d[3] = {0,0,0};
	
	GetCouplingMatrix(M, d, f, index.data(), S, flipping_S, pb, grid);

	gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);

	double det = abs( LU[0][0] * LU[1][1] * LU[2][2]);
	free_matrix(M,grid.dim-1,grid.dim-1);
	free_matrix(LU,grid.dim-1,grid.dim-1);
	free_matrix(f,grid.dim-1,4,4,4);

	return det;
}
// interpolate to get f at ghost point
void SetF(Index ind, double ***b, double*** S, double*** sign_surf, PBData& pb, GridData &grid){

	if (IsGhost(ind,S,sign_surf)){
		bool found = false;	
		for (int r = 0; r < 3; r++){
			for( int s = -1; s <= 1; s+=2){
				Index nbr1 = UnitIncrement(ind,r,s);
				Index nbr2 = UnitIncrement(ind,r,s);
				if(!IsGhost(nbr1,S,sign_surf) && !IsGhost(nbr2,S,sign_surf) 
					&& SameSide(sign_surf, nbr1, ind) && SameSide(sign_surf, nbr2, ind)){
					// two nbr not ghost and same side
					double f1 = getf(nbr1.data(),0,0,0.0,GetSign(sign_surf,nbr1),pb,grid);
					double f2 = getf(nbr1.data(),0,0,0.0,GetSign(sign_surf,nbr2),pb,grid);
					setvalarray(b,ind.data(),(f1+f2)/2.0);
					found = true;
					return;
				}
			}
		}
		if(found==false){
			fprintf(stderr, "Fail to interpolate f at (%i,%i,%i)\n", ind[0],ind[1],ind[2]);
			exit(1);
		}
	}else{
		// if not ghost
		if (evalarray(sign_surf,ind.data()) < 0.0)
			setvalarray(b,ind.data(),getf(ind.data(),0,0,0.0,-1,pb,grid));
		else
			setvalarray(b,ind.data(),getf(ind.data(),0,0,0.0,1,pb,grid));
	}

	return;
}

// lineary system for icim
void linearsystem_icim(SparseElt2**** &A, double ***b, double ***S, double*** sign_surf, PBData &pb, GridData &grid)
{

	int count[4] = {0};

	cout << "in linearsystem_icim" << endl;

	clearsparse(A,grid);

	// init b
	for(int i = 0 ; i <= grid.nx[0]; i++){
		for(int j = 0 ; j <= grid.nx[1]; j++){
			for(int k = 0 ; k <= grid.nx[2]; k++){
				Index tindex = {i,j,k};
				// SetF(tindex,b,S,sign_surf,pb, grid);
				setvalarray(b,tindex,getf(tindex.data(),0,0,0.0, GetSign(sign_surf, tindex),pb,grid));
			}
		}
	}

	// finite different at each grid point
	for(int i = 0 ; i <= grid.nx[0]; i++){
		for(int j = 0 ; j <= grid.nx[1]; j++){
			for(int k = 0 ; k <= grid.nx[2]; k++){
				Index tindex = {i,j,k};
				char thestatus = GetType( tindex, sign_surf, grid );
				
				Index rindex = tindex;


				if (thestatus == 0){
					// interior points
					interiorptsmall(A, b, tindex.data(), S, pb, grid);
					(count[1])++;
				
				}else if (thestatus == 4){
					// boundary points
					sparse2(tindex.data(), tindex.data(), A, 1.0, grid);
					
					for (int r = 0; r < grid.dim; r++){
						array<double,3> x = sub2coord(tindex, grid);
						setvalarray(b, tindex.data(), DBC(x.data(), grid.dim, 0.0) );
					}
					(count[0])++;
				
				}else{

					icim2all(A, b, tindex.data(), S, sign_surf, pb, grid);
					// cim2(A, b, Dusmall, buildsize, tindex.data(), gamma, S, pb, grid);

					(count[2])++;
				}
			}
		}
	}


	cout << "found " << count[0] << " boundary pts, " << count[1] << " interior pts, " << count[2] << " cim pts" <<endl;
}

// has two nbr not ghost
bool Exist2NotGhostNbr(int &dim, int &dir, Index ind, double *** S, double *** sign_surf, GridData& grid){
	for(int r = 0; r < 3; r++){
		for(int s = -1; s<=1; s+=2){
			Index step1 = UnitIncrement(ind, r, s);
			Index step2 = UnitIncrement(ind, r, s*2);
			if(!outofbound( step2, grid)){
				if (  SameSide(S, step1, ind) && (!IsGhost(step1, S, sign_surf)) &&
					  SameSide(S, step2, ind) && (!IsGhost(step2, S, sign_surf))    ) {
					dim = r;
					dir = s;
					return true;
				}
			}else{
				continue;
			}
		}
	}
	return false;
}

// find nearest nbr not that is not ghost amd same side
Index GetNearestRealNbr(Index index, double *** sign_surf, double *** S,  GridData &grid){
	
	Index nearest;

	vector<Index> nbr_order = GetCubeNeighbors(index, grid, CUBENBR);
	nbr_order.insert(nbr_order.begin(),index);

	std::sort( nbr_order.begin(), nbr_order.end(),
    [index](Index &a,  Index &b){ return ( Norm2(Minus(a,index)) < Norm2(Minus(b,index)) );});  

	for(auto temp : nbr_order){
		if (!IsGhost(temp,sign_surf,S) && SameSide(S,temp,index)){
			nearest = temp;
			break;
		}
	}
	return nearest;

}

// find nearest nbr not that is not ghost amd same side and is interface point
Index GetNearestGhostInterfaceNbr(Index index, double *** S,  GridData &grid){
	
	Index nearest;
	vector<Index> nbr_order = GetCubeNeighbors(index, grid, CUBENBR);
	nbr_order.insert(nbr_order.begin(),index);

	std::sort( nbr_order.begin(), nbr_order.end(),
    [index](Index &a,  Index &b){ return ( Norm2(Minus(a,index)) < Norm2(Minus(b,index)) );});  

	for(auto temp : nbr_order){

		int ind1d = sub2ind(temp.data(), grid.nx, grid.dim);

		if ( !(D2uMap.find(ind1d)==D2uMap.end()) && D2uMap[ind1d].sign * evalarray(S, index) > 0.0 ){
			// exist in map and same sign
			nearest = temp;
			break;
		}
	}
	return nearest;

}


void RecoverD2u(Index index, double *** u, GridData &grid){

	int ind1d = sub2ind(index.data(),grid.nx,grid.dim);

	if(D2uMap.find( ind1d  )==D2uMap.end()){
		cerr<<"index not in map"<<endl;
		exit(1);
	}

	Info info = D2uMap[ind1d];

	array<double,3> uxx;
	int N = 4;


	double d[3];
	double** M = matrix(grid.dim-1,grid.dim-1);
	double**** f = matrix(grid.dim-1,N,N,N);

	GetInfo(M,d,f,info);

	// if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
	// 	printf("[RecoverD2u] (%i,%i,%i)\n", index[0],index[1],index[2]);		
	// 	double v[3] = {0};
	// 	double tempuxcoef[3] = {0};
	// 	double tempuxxcoef[3] = {0};
		
	// 	cout<<"recover M matrix"<<endl;
	// 	printMat(3,3,M);
	// 	cout<<"recover d vector"<<endl;
	// 	printMat(3,d);
	// }

	int PLR[grid.dim], PLC[grid.dim];
	double** LU = matrix(grid.dim-1,grid.dim-1);

	gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);
	
	forwardbacksub0(d,d,LU,PLR,PLC,grid.dim-1);


	// copy from checkDanswer2
	for (int r = 0; r < grid.dim; r++)
		uxx[r] = d[r]/(grid.dx[r]*grid.dx[r]);
	 
	int tindex[3],rindex[3];
	double temp[3];
	for (int i = 0; i < grid.dim; i++)
		tindex[i] = 0;
	while (tindex[0] < 5) //index is local index for u
	{
		for (int r = 0; r < grid.dim; r++)
			temp[r] = evalarray(f[r],tindex);//temp is the 3x1 vector, one column of D, corresponding to u_tindex
		forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
		for (int r = 0; r < grid.dim; r++)
			rindex[r] = index[r]-2+tindex[r];//rindex is global index for u
		for (int r = 0; r < grid.dim; r++)
			uxx[r] += temp[r]/(grid.dx[r]*grid.dx[r])*evalarray(u,rindex);//calculate inv(M)

		(tindex[grid.dim-1])++;
		for (int i = grid.dim-1; i > 0 && tindex[i] >= 5; i--)
		{
			tindex[i] = 0;
			(tindex[i-1])++;
		}
	}


	free_matrix(M,grid.dim-1,grid.dim-1);
	free_matrix(LU,grid.dim-1,grid.dim-1);
	free_matrix(f,grid.dim-1,4,4,4);
	
	D2uMap[ind1d].d2u[0] = uxx[0];
	D2uMap[ind1d].d2u[1] = uxx[1];
	D2uMap[ind1d].d2u[2] = uxx[2];

}


double CentralDiffD2u(Index ind, double*** u, int dim, GridData &grid){
	return 1.0 / pow(grid.dx[dim],2) * (
		evalarray(u, UnitIncrement(ind, dim, 1))
		- 2.0 * evalarray(u, ind)
		+ evalarray(u, UnitIncrement(ind, dim, -1)));
}

// find D2u by looking up nearest interface in sign_surf
array<double,3> ComputeD2uLookup(Index index, double*** S, GridData& grid){
	
	Vector3d D2u = {NAN,NAN,NAN};

	Index ind_to_eval = GetNearestGhostInterfaceNbr(index, S, grid);
	D2u[0] = D2uMap[sub2ind(ind_to_eval.data(), grid.nx, grid.dim)].d2u[0];
	D2u[1] = D2uMap[sub2ind(ind_to_eval.data(), grid.nx, grid.dim)].d2u[1];
	D2u[2] = D2uMap[sub2ind(ind_to_eval.data(), grid.nx, grid.dim)].d2u[2];
	
	if (ind_to_eval.empty()){
		fprintf(stderr, "(%i, %i, %i) Can not find same side non-ghost nbr\n", index[0], index[1], index[2]);
	}

	// error of uxx
	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		printf("[ComputeD2uLookup] (%i,%i,%i), eval at (%i,%i,%i)\n",index[0],index[1],index[2],ind_to_eval[0],ind_to_eval[1],ind_to_eval[2]);
		double exact_uxx[3];
		for(int j = 0; j < 3; j++){
			exact_uxx[j] = getD2u(index.data(),j,j,0,0,0,GetSign(S,index),grid);
			printf("D2u[%i][%i] compute = %f, exact = %f, err = %f\n", j, j, D2u[j], exact_uxx[j], abs(D2u[j] - exact_uxx[j]) );
		}
	}
		
	
	return D2u;
}


array<double,3> ComputeD2u(Index index, double***u, double*** S, GridData& grid){
	Vector3d D2u = {NAN,NAN,NAN};
	Index g = Indicator(index, S, grid);

	Index ind_to_eval = GetNearestGhostInterfaceNbr(index, S, grid);

	if(index == ind_to_eval){
		D2u[0] = D2uMap[sub2ind(ind_to_eval.data(), grid.nx, grid.dim)].d2u[0];
		D2u[1] = D2uMap[sub2ind(ind_to_eval.data(), grid.nx, grid.dim)].d2u[1];
		D2u[2] = D2uMap[sub2ind(ind_to_eval.data(), grid.nx, grid.dim)].d2u[2];
	}else{
		for(int r : {0, 1, 2}){
			if(g[r] == 2){
				D2u[r] = CentralDiffD2u(index, u, r, grid);
			}else{
				D2u[r] = D2uMap[sub2ind(ind_to_eval.data(), grid.nx, grid.dim)].d2u[r];
			}
		}	
	}

	// error of uxx
	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		printf("[ComputeD2uLookup] (%i,%i,%i), eval at (%i,%i,%i)\n",index[0],index[1],index[2],ind_to_eval[0],ind_to_eval[1],ind_to_eval[2]);
		double exact_uxx[3];
		for(int j = 0; j < 3; j++){
			exact_uxx[j] = getD2u(index.data(),j,j,0,0,0,GetSign(S,index),grid);
			printf("D2u[%i][%i] compute = %f, exact = %f, err = %f\n", j, j, D2u[j], exact_uxx[j], abs(D2u[j] - exact_uxx[j]) );
		}
	}
		
	
	return D2u;
}




void ComputeMixD2u(double** D2u, double*** u_ghost, double *** sign_surf, double *** S, Index index, GridData &grid){
	int mid = 2;
	int N = 2 * mid;
	double *** mixDu_ucoef = matrix(N, N, N); // mixDu in r, rstar plane in terms of u
	double mixDu_uxxcoef[3] = {0}; // mixDu in r, rstar plane in terms of uxx

	array<double,3> uxx = ComputeD2uLookup( index, S, grid);

	Index local_mid = {mid,mid,mid};
	for(int m = 0; m < 3; m ++){
		for(int n = m + 1; n < 3; n ++){
			D2u[m][n] = 0.0;
			if (HasMixD2u(mixDu_ucoef,  mixDu_uxxcoef, m, n, index , sign_surf, grid)){
				//eval 
				for( int i = 0; i <= N; i ++){
					for( int j = 0; j <= N; j ++){
						for( int k = 0; k <= N; k ++){
							Index tindex = {i,j,k};
							Index globindex = Add(index, Minus(tindex, local_mid));
							D2u[m][n] += evalarray(u_ghost,globindex) * mixDu_ucoef[i][j][k];
						}
					}
				}

				for(int i = 0; i < 3; i++){
					D2u[m][n] += mixDu_uxxcoef[i] * uxx[i];
				}


				if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
					printf("[ComputeMixD2u] (%i,%i,%i)\n",index[0],index[1],index[2]);
					double exact_uxx = getD2u(index.data(),m,n,0,0,0,GetSign(S,index),grid);
					printf("D2u[%i][%i] compute = %f, exact = %f, err = %f\n", m, n, D2u[m][n], exact_uxx, abs(D2u[m][n] - exact_uxx));

					printf("u coeff\n");
					PrintCoef(mixDu_ucoef,N);
					printf("uxx coeff\n");
					PrintCoef(mixDu_uxxcoef,2);
				}
				

			}else{
				fprintf(stderr, "cannot approximate mixDu at plane (%i,%i) at (%i, %i, %i)\n", m, n, index[0], index[1], index[2]);
				exit(1);
			}		
		}
	}

	for(int m = 0; m < 3; m ++){
		for(int n = m + 1; n < 3; n ++){
			D2u[n][m] = D2u[m][n];
		}
	}


	free_matrix(mixDu_ucoef,N,N,N);
	return ;
}


void ComputeDuGridpoint(double Du[], double*** u_ghost, double *** sign_surf, double *** S, Index index, GridData &grid){
	int mid = 2;
	int N = 2 * mid;
	Index local_mid = {mid,mid,mid};


	double**** D1ucoef = new double ***[grid.dim]; 
	for (int r = 0; r < grid.dim; r++)
		D1ucoef[r] = matrix(N,N,N);
	double** D1uxxcoef = matrix(grid.dim-1,grid.dim-1);


	GetIcimDuGridpoint( D1ucoef, D1uxxcoef, sign_surf, index.data(), grid);

	// get D2u get grid point
	array<double,3> uxx = ComputeD2u( index, u_ghost, S, grid);

	for(int r = 0; r < 3; r ++){
		// evalute coeff with computed data
		Du[r] = 0;
		for( int i = 0; i <= N; i ++){
			for( int j = 0; j <= N; j ++){
				for( int k = 0; k <= N; k ++){
					Index tindex = {i,j,k};
					Index globindex = Add(index, Minus(tindex, local_mid));

					Du[r] += evalarray(u_ghost, globindex) * D1ucoef[r][i][j][k];
				}
			}
		}

		for(int j = 0; j < 3; j ++){
			Du[r] += D1uxxcoef[r][j] * uxx[j];
		}
		
	}

		// error of uxx
	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		printf("[ComputeDuGridpoint] (%i,%i,%i)\n",index[0],index[1],index[2]);
		double exact_du[3];
		for(int r = 0; r < 3; r ++){
			exact_du[r] = getDu(index.data(), r, 0, 0, 0, GetSign(S, index), grid);
			printf("Du[%i] compute = %f, exact = %f, err = %f\n", r, Du[r], exact_du[r], abs(Du[r] - exact_du[r]) );
			
		}
	}

	return;

}

double TaylorExpand(double u0, double* ux, double** uxx, double* h){
	double uval = u0;

	if(ux){
		for(int r = 0; r < 3; r ++){
			uval += h[r] * ux[r];
		}	
	}
	
	if(uxx){
		for(int m = 0; m < 3; m ++){
			for(int n = 0; n < 3; n ++){
				uval +=  1.0/2.0 * h[m] * h[n] * uxx[min(m,n)][max(m,n)];	
			}
		}	
	}
	
	return uval;
}


// simple way to reconstruct u at index
// use central diff for D2u, need two consecutive non-ghost nbr
// icim paper eqn 18
double ComputeSimple(double*** u_ghost, double *** sign_surf, double *** S, Index index, GridData &grid){
	for(int r : {0, 1, 2}){
		for(int s : {-1, 1}){
			Index step1 = UnitIncrement(index, r, s);
			Index step2 = UnitIncrement(index, r, 2*s);
			if(SameSide(S, step1, index) && SameSide(S, step2, index) && !IsGhost(step1, S, sign_surf) && !IsGhost(step2, S, sign_surf) ){
				array<double,3> uxx = ComputeD2uLookup( step1,  S, grid);
				return 2 * evalarray(u_ghost,step1) - evalarray(u_ghost,step2) + pow(grid.dx[r], 2) * uxx[r];
			}
		}
	}
	return NAN;
}



// compute u and du value at interface or grid point by extrapolation from nearby non-ghost points
void ComputeByExtrapolate(double& uval, Vector3d& du, 
 double*** u_ghost, double *** sign_surf, double *** S, Index index, int rstar, int sstar, double alpha, GridData &grid){
	int mid = 2;
	int N = 2*mid;
	Index local_mid = {mid, mid, mid};
	
	// add nbr, sort by distance to index
	Index nearest = GetNearestRealNbr(index, sign_surf, S, grid);

	if (nearest.empty()){
		cerr<<"no neighbor found"<<endl;
		return;
	}

	// compute u, du, d2u at nearest grid point
	double u_grid;
	double Du_grid[3];
	double** D2u = matrix(2,2);

	u_grid = evalarray(u_ghost, nearest);
	// second deri at nearest
	array<double,3> uxx = ComputeD2u(nearest, u_ghost, S, grid);

	ComputeMixD2u(D2u, u_ghost, sign_surf, S, nearest, grid);

	// fill in diagonal
	D2u[0][0] = uxx[0];
	D2u[1][1] = uxx[1];
	D2u[2][2] = uxx[2];

	ComputeDuGridpoint(Du_grid, u_ghost, sign_surf, S, nearest, grid);	


	// step from nearest to index
	Index s = Minus(index, nearest);
	// h is distance from nearest to interface point
	Vector3d h = {s[0] * grid.dx[0], s[1] * grid.dx[1], s[2] * grid.dx[2]};

	// get distance to inteface
	h[rstar] += sstar * alpha * grid.dx[rstar];

	for(int i = 0; i < 3; i ++){
		du[i] = TaylorExpand(Du_grid[i], D2u[i], nullptr, h.data());
	}
	// taylor series
	uval = TaylorExpand(u_grid, Du_grid, D2u, h.data());
	

	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){

		printf("[ComputeByExtrapolate] (%i,%i,%i), rstar = %i, sstar = %i, alpha = %f, nearest = (%i,%i,%i), h = (%f,%f,%f)\n",
			index[0],index[1],index[2],rstar, sstar, alpha, nearest[0],nearest[1],nearest[2], h[0],h[1],h[2]);
		
		
		// check nearest u error
		double nearest_exact_u = getu(nearest.data(), 0, 0, 0, GetSign(S, nearest), grid);
		PrintError("nearest u exact", nearest_exact_u, "nearest u compute",evalarray(u_ghost,nearest), "err");
		
		// check nearest du error
		double nearest_exact_du[3];
		for(int r = 0; r < 3; r ++){
			nearest_exact_du[r] = getDu(nearest.data(), r, 0, 0, 0, GetSign(S, nearest), grid);
		}		
		PrintError(3, "nearest du exact", nearest_exact_du, "nearest du compute", Du_grid, "loc truc err");

		// check nearest d2u error 
		double** nearest_exact_d2u = matrix(2,2);
		for(int m = 0; m < 3; m++){
			for(int n = 0; n < 3; n++){
				nearest_exact_d2u[m][n] = getD2u(nearest.data(), m, n, 0, 0, 0, GetSign(S, nearest), grid);
			}
			PrintError(3, "nearest d2u exact", nearest_exact_d2u[m], "nearest d2u compute", D2u[m], "err");
		}

		// check Du error
		double exact_du[3];
		for(int r = 0; r < 3; r ++){
			exact_du[r] = getDu(index.data(), r, rstar, sstar, alpha, GetSign(S, index), grid);
		}
		double Du_taylor[3];
		for(int i = 0; i < 3; i ++){
			Du_taylor[i] = TaylorExpand(nearest_exact_du[i], nearest_exact_d2u[i], nullptr, h.data());
		}

		PrintError(3, "du exact", exact_du, "du taylor with exact", Du_taylor, "loc truc err");
		PrintError(3, "du exact", exact_du, "du compute", du.data(), "err");

		// check u error
		double exact_u = getu(index.data(),rstar,sstar,alpha,GetSign(S,index),grid);
		double uval_taylor = TaylorExpand(nearest_exact_u, nearest_exact_du, nearest_exact_d2u, h.data());
		PrintError("u exact", exact_u, "u taylor with exact", uval_taylor, "loc trunc err");
		PrintError("u exact", exact_u, "u compute", uval, "err");
		cout<<endl;

		free_matrix(nearest_exact_d2u,2,2);
		
	}

	free_matrix(D2u,2,2);
}

// check du at interface
void CheckIcimDu(double*** u_ghost, double *** u_real, double*** sign_surf, double ***S,  PBData &pb, GridData &grid){
	cout<<"[Check icim Du]"<<endl;
	double max_err = 0.0;
	Index max_idx;
	int nancount = 0; // number of interface points fail
	vector<double> err_all;


	for(int i = 0; i <= grid.nx[0]; i ++){
	for(int j = 0; j <= grid.nx[1]; j ++){
	for(int k = 0; k <= grid.nx[2]; k ++){
		Index index = {i, j, k};
		if(nearinterface(index, S, grid)){
			// For interface points
			for(int rstar = 0; rstar < 3; rstar ++){
				for(int sstar = -1; sstar <=1; sstar+=2){
					// get direction of interface
					Index rindex = UnitIncrement(index, rstar, sstar);
					if(!SameSide(S, index, rindex )){
						double alpha;
						double tangent[3], normal[3];
						Vector3d exact_Du = {0,0,0};
						
						getinterfaceinfo(alpha,tangent,normal,S,index.data(),rindex.data(),grid);

						// get excat Du at itnerface
						for(int r = 0; r < 3; r ++){
							exact_Du[r] = getDu(index.data(), r, rstar, sstar, alpha, GetSign(S, index), grid);
						}

						double dummy;
						Vector3d apprx_Du;
						ComputeByExtrapolate(dummy, apprx_Du, u_real, S, S, index,  rstar,  sstar, alpha, grid);

						if( isnan(apprx_Du)){
							nancount++;
							continue;
						}

						for(int r = 0; r < 3; r ++){
							err_all.push_back(abs(apprx_Du[r] - exact_Du[r]));
						}

						double err = NormInf(Minus(apprx_Du, exact_Du));

						if(err > max_err){
							max_err = err;
							max_idx = index;
						}
													

					}	
				}
			}
		}

	}
	}
	}
	
	printf("fail to extrapolate = %i \n", nancount);
	printf("Du max err = %.8f at %i %i %i\n", max_err, max_idx[0], max_idx[1], max_idx[2]);
	printf("Du rmse = %.8f \n", rmse(err_all));

}


// check u at interface
void CheckIcimU(double*** u_ghost, double *** u, double*** sign_surf, double ***S,  PBData &pb, GridData &grid){
	cout<<"[Check icim u at interface]"<<endl;
	double max_err = 0.0;
	Index max_idx;

	vector<double> err_all;


	for(int i = 0; i <= grid.nx[0]; i ++){
	for(int j = 0; j <= grid.nx[0]; j ++){
	for(int k = 0; k <= grid.nx[0]; k ++){
		Index index = {i, j, k};
		if(nearinterface(index, S, grid)){
			// For interface points
			for(int rstar = 0; rstar < 3; rstar ++){
				for(int sstar = -1; sstar <=1; sstar+=2){
					// get direction of interface
					Index rindex = UnitIncrement(index, rstar, sstar);
					if(!SameSide(S, index, rindex )){
						double alpha;
						double tangent[3], normal[3];
						
						getinterfaceinfo(alpha,tangent,normal,S,index.data(),rindex.data(),grid);

						// get excat u at itnerface
						double exact_u = getu(index.data(), rstar, sstar, alpha, GetSign(S, index), grid);
						
						double uval;
						Vector3d dummy;
						ComputeByExtrapolate(uval, dummy, u_ghost, sign_surf, S, index,  rstar,  sstar, alpha, grid);

						double err = abs(uval - exact_u);

						err_all.push_back(err);

						if(err > max_err){
							max_err = err;
							max_idx = index;
						}
													

					}	
				}
			}
		}

	}
	}
	}

	printf("u max err = %.8f at %i %i %i\n", max_err, max_idx[0], max_idx[1], max_idx[2]);
	printf("u rmse = %.8f \n", rmse(err_all));

}


void CheckErrGrid(double *** u, double ***S,  PBData &pb, GridData &grid){
	cout<<"[Check icim u grid point]"<<endl;
	double max_err = 0.0;
	Index max_idx;

	vector<double> err_all;
	int nancount = 0;

	for(int i = 0; i <= grid.nx[0]; i ++){
	for(int j = 0; j <= grid.nx[0]; j ++){
	for(int k = 0; k <= grid.nx[0]; k ++){
		Index index = {i, j, k};
		if(isnan(evalarray(u,index))){
			nancount++;
			continue;
		}

		double err = abs( evalarray(u,index)-getu(index.data(),0,0,0.0,GetSign(S,index),grid));

		err_all.push_back(err);

		if(err > max_err){
			max_err = err;
			max_idx = index;
		}

	}
	}
	}

	printf("u max err = %.8f at %i %i %i\n", max_err, max_idx[0], max_idx[1], max_idx[2]);
	printf("u rmse = %.8f \n", rmse(err_all));
	printf("fail to recontruct u nan = %i \n", nancount);
	
}


// get u and ux on the other side of interface
void GetRowRhs(double &u, double &h, Vector3d& ux, Vector3d& normal, double& sigma, double& tau,
 double*** u_ghost, Index index, double *** sign_surf, double *** S, int r, int s, PBData &pb, GridData& grid){

	int step = s;
	Vector3d xindex = sub2coord(index, grid);
	while(true){
		Index idx_other = UnitIncrement(index, r, step);
		if(IsGhost(idx_other, sign_surf, S)){
			step += s;
			continue;
		}else{
			// find a non-ghost point
			if(SameSide(S, index, idx_other)){
				u = evalarray(u_ghost, idx_other);
				h = (idx_other[r] - index[r]) * grid.dx[r];
				
				sigma = NAN;
				tau = NAN;
				ux = Vector3d {3,NAN};


			}else{
				double alpha;
				Vector3d tangent;

				Index idx_same = UnitIncrement(idx_other, r, -s);//step back

				if(SameSide(S, idx_same, idx_other)){
					cerr<<"no interface between this sie and  idx_other"<<endl;
					exit(1);
				}
				getinterfaceinfo(alpha, tangent.data(), normal.data(), S, idx_same.data(), idx_other.data(), grid);

				Vector3d x = sub2coord(idx_same, grid);
				x[r] += s * grid.dx[r] * alpha;// interface location
				h = x[r] - xindex[r]; // distance to interface

				sigma = getsigma(x.data(), normal.data(), pb, grid);
				tau = gettau(x.data(), pb, grid);

				ComputeByExtrapolate(u, ux, u_ghost, sign_surf, S, idx_other, r, -s, 1.0 - alpha, grid);
			}
			break;// break out of while loop
		}

		if(atbound(index, grid)){
			fprintf(stderr, "get to boundary at (%i,%i,%i) r = %i, s = %i \n", index[0], index[1], index[2], r, s);
			exit(1);
		}
	}

}


int D2uCol(int m, int n){
	if(m == n){
		return 4 + m;
	}else{
		// mixded du
		return 6 + min(m,n) + max(m,n);
	}
}

int DuCol(int m){
	return m + 1;
}



void ComputeIcimLinSolve(double& u, double* du, double** d2u,
 double*** u_ghost, double *** sign_surf, double *** S, Index index, PBData &pb, GridData &grid){
	const int rownum = 10;
 	// for debug
 	double xexact[rownum];
	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		// check each row
		xexact[0] = getu(index.data(), 0, 0, 0, GetSign(S, index), grid);
		// check nearest du error
		for(int r = 0; r < 3; r ++){
			xexact[DuCol(r)] = getDu(index.data(), r, 0, 0, 0, GetSign(S, index), grid);
		}		
		// check nearest d2u error 
		for(int m = 0; m < 3; m++){
			for(int n = m; n < 3; n++){
				xexact[D2uCol(m,n)] = getD2u(index.data(), m, n, 0, 0, 0, GetSign(S, index), grid);
			}
		}
		PrintFlipSurf(sign_surf, S, index.data(), 2);
	}


	
	Vector3d grid_coord = sub2coord(index, grid);
	vector< vector<double> > A( rownum, vector<double>(rownum, 0.0));
	vector<double> b(rownum, 0.0);

	// first row , from pde
	int row = 0;
	double thesign = GetSign(S, index);
	double ehere = (thesign < 0)? pb.epsilonm : pb.epsilonp;
	double ethere = (thesign < 0)? pb.epsilonp : pb.epsilonm;
	A[0][4] = - ehere;
	A[0][5] = - ehere;
	A[0][6] = - ehere;
	b[0] = getf(grid_coord.data(), GetSign(S, index), pb, grid);

	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		// check each row
		double lhs = getdotprod(xexact, (A[row]).data(), 10);
		printf("row %i, lhs = %f, rhs = %f, err = %f\n", row, lhs, b[row], abs(lhs - b[row]));
	}


	row++;

	// sort dim by descending number of interface
	Index g = Indicator(index, S, grid);
	Index dim_order = {0,1,2};
	sort(dim_order.begin(),dim_order.end(), [g](int i1, int i2) { return g[i1] < g[i2]; } );
	list<pair<int,int>> r_s_ordered;
	for(int k = 0; k < 3; k++){
		r_s_ordered.push_back(make_pair(dim_order[k],-1));
		r_s_ordered.push_back(make_pair(dim_order[k],1));
	}


	while(!r_s_ordered.empty()){
		int r = r_s_ordered.front().first;
		int s = r_s_ordered.front().second;
		r_s_ordered.pop_front();
		
		Vector3d ux,  normal;

		double u, h, sigma, tau;
		GetRowRhs(u, h, ux, normal, sigma, tau,  u_ghost, index, sign_surf, S, r, s, pb, grid);		
		if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		printf("normal = (%f,%f,%f), sigma = %f, tau =  %f \n",normal[0], normal[1], normal[2], sigma, tau);
		}

		bool sameside = isnan(sigma);

		if(sameside && !isnan(u)){
			Index othergrid = UnitIncrement(index, r, s);
			A[row][0] = 1;
			A[row][DuCol(r)] = h;
			A[row][D2uCol(r,r)] = 0.5 * h * h;
			b[row] = u;

			if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
				double lhs = getdotprod(xexact, (A[row]).data(), 10);
				double exactrhs = getu(othergrid.data(), r, s, 1.0, GetSign(S, othergrid.data()), grid);
				printf("row %i, r = %i, s = %i, u, h = %f, lhs = %f, rhs compute= %f, residual = %f, rhs exact = %f , rhs error = %f \n",
				 row, r, s, h, lhs, b[row], abs(lhs - b[row]), exactrhs, abs(b[row] - exactrhs));
			}

			row++;
			if(row == rownum){
				break;
			}
		}


		if(!sameside){
			Vector3d interface_coord = grid_coord;
			interface_coord[r] += h;

			if(!isnan(u)){
				// u equation
				A[row][0] = 1;
				A[row][DuCol(r)] = h;
				A[row][D2uCol(r,r)] = 0.5 * h * h;
				b[row] = u + thesign * tau;

				if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){

					double lhs = getdotprod(xexact, (A[row]).data(), 10);
					double exactrhs = getu(interface_coord.data(),GetSign(S, index.data()), grid);
					printf("row %i, r = %i, s = %i, u, h = %f, lhs = %f, rhs compute= %f, residual = %f \n",
					 row, r, s, h, lhs, b[row], abs(lhs - b[row]));
					// printf("computed u other side = %f, tau = %f, rhs exact = %f , rhs error = %f \n",
					//  u, tau, exactrhs, abs(b[row] - exactrhs));
				}


				row++;
				if(row == rownum){
					break;
				}
			}

			if(all_of(ux.cbegin(), ux.cend(), [](double e) { return !isnan(e); })){
				// du equation
				A[row][DuCol(0)] = ehere * normal[0]; // ux
				A[row][DuCol(1)] = ehere * normal[1];	// uy
				A[row][DuCol(2)] = ehere * normal[2];	// uz

				A[row][D2uCol(r,0)] = ehere * normal[0] * h; // ux
				A[row][D2uCol(r,1)] = ehere * normal[1] * h; // uy
				A[row][D2uCol(r,2)] = ehere * normal[2] * h; // uz

				b[row] = ethere * (Dot(normal, ux)) + thesign * sigma;

				if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
					
					double lhs = getdotprod(xexact, (A[row]).data(), 10);
					double exactrhs = 0.0; // epsilon- dot(grad(u-),n)
					for(int r = 0; r < 3; r ++){
						exactrhs += ehere * normal[r] * getDu( interface_coord.data(), r, GetSign(S, index), grid);
					}		
					printf("row %i, r = %i, s = %i, Du, h = %f, lhs = %f, rhs compute= %f, residual = %f, rhs exact = %f , rhs error = %f \n",
					 row, r, s, h, lhs, b[row], abs(lhs - b[row]), exactrhs, abs(b[row] - exactrhs));
					// printf("normal = (%f,%f,%f), sigma = %f, computed other side Dot(normal, ux) = %f \n",
					//  normal[0], normal[1], normal[2], sigma, Dot(normal, ux));
				}


				row++;
				if(row == rownum){
					break;
				}
			}
		}
	}

	

	if(row < 10){
		fprintf(stderr, "Not enougth equations at (%i,%i,%i)\n", index[0],index[1],index[2]);
		SetMatrix(d2u,(double)NAN, 2, 2);
		SetMatrix(du, (double)NAN, 2);
	}

	vector<double> x = gesolve(A, b);
	u = x[0];

	for( int m = 0; m < 3; m ++){
		du[m] = x[DuCol(m)];
		for(int n = m; n < 3; n++){
				d2u[m][n] = x[D2uCol(m,n)];	
		}
	}

	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		printMat(10,10,A);
		printMat(10,b);
		PrintError(10, "solved", x.data(), "exact", xexact, "err");
	}
	
}



void ComputeIcimLinSolve7var(double& u, double* du, double** d2u,
 double*** u_ghost, double *** sign_surf, double *** S, Index index, PBData &pb, GridData &grid){
	const int rownum = 7;
 	// for debug
 	double xexact[rownum];
	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		// check each row
		xexact[0] = getu(index.data(), 0, 0, 0, GetSign(S, index), grid);
		// check nearest du error
		for(int r = 0; r < 3; r ++){
			xexact[DuCol(r)] = getDu(index.data(), r, 0, 0, 0, GetSign(S, index), grid);
		}		
		// check nearest d2u error 
		for(int m = 0; m < 3; m++){
			for(int n = m; n < 3; n++){
				xexact[D2uCol(m,n)] = getD2u(index.data(), m, n, 0, 0, 0, GetSign(S, index), grid);
			}
		}
		PrintFlipSurf(sign_surf, S, index.data(), 2);
	}


	
	Vector3d grid_coord = sub2coord(index, grid);
	vector< vector<double> > A( rownum, vector<double>(rownum, 0.0));
	vector<double> b(rownum, 0.0);

	// first row , from pde
	int row = 0;
	double thesign = GetSign(S, index);
	double ehere = (thesign < 0)? pb.epsilonm : pb.epsilonp;
	double ethere = (thesign < 0)? pb.epsilonp : pb.epsilonm;
	A[0][4] = - ehere;
	A[0][5] = - ehere;
	A[0][6] = - ehere;
	b[0] = getf(grid_coord.data(), GetSign(S, index), pb, grid);

	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		// check each row
		double lhs = getdotprod(xexact, (A[row]).data(), 10);
		printf("row %i, lhs = %f, rhs = %f, err = %f\n", row, lhs, b[row], abs(lhs - b[row]));
	}


	row++;

	// sort dim by descending number of interface
	Index g = Indicator(index, S, grid);
	Index dim_order = {0,1,2};
	sort(dim_order.begin(),dim_order.end(), [g](int i1, int i2) { return g[i1] < g[i2]; } );
	list<pair<int,int>> r_s_ordered;
	for(int k = 0; k < 3; k++){
		r_s_ordered.push_back(make_pair(dim_order[k],-1));
		r_s_ordered.push_back(make_pair(dim_order[k],1));
	}


	while(!r_s_ordered.empty()){
		int r = r_s_ordered.front().first;
		int s = r_s_ordered.front().second;
		r_s_ordered.pop_front();
		
		Vector3d ux,  normal;

		double u, h, sigma, tau;
		GetRowRhs(u, h, ux, normal, sigma, tau,  u_ghost, index, sign_surf, S, r, s, pb, grid);		
		if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		printf("normal = (%f,%f,%f), sigma = %f, tau =  %f \n",normal[0], normal[1], normal[2], sigma, tau);
		}

		bool sameside = isnan(sigma);

		if(sameside && !isnan(u)){
			Index othergrid = UnitIncrement(index, r, s);
			A[row][0] = 1;
			A[row][DuCol(r)] = h;
			A[row][D2uCol(r,r)] = 0.5 * h * h;
			b[row] = u;

			if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
				double lhs = getdotprod(xexact, (A[row]).data(), 10);
				double exactrhs = getu(othergrid.data(), r, s, 1.0, GetSign(S, othergrid.data()), grid);
				printf("row %i, r = %i, s = %i, u, h = %f, lhs = %f, rhs compute= %f, residual = %f, rhs exact = %f , rhs error = %f \n",
				 row, r, s, h, lhs, b[row], abs(lhs - b[row]), exactrhs, abs(b[row] - exactrhs));
			}

			row++;
			if(row == rownum){
				break;
			}
		}


		if(!sameside){
			Vector3d interface_coord = grid_coord;
			interface_coord[r] += h;

			if(!isnan(u)){
				// u equation
				A[row][0] = 1;
				A[row][DuCol(r)] = h;
				A[row][D2uCol(r,r)] = 0.5 * h * h;
				b[row] = u + thesign * tau;

				if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){

					double lhs = getdotprod(xexact, (A[row]).data(), 10);
					double exactrhs = getu(interface_coord.data(),GetSign(S, index.data()), grid);
					printf("row %i, r = %i, s = %i, u, h = %f, lhs = %f, rhs compute= %f, residual = %f \n", row, r, s, h, lhs, b[row], abs(lhs - b[row]));
					// printf("computed u other side = %f, tau = %f, rhs exact = %f , rhs error = %f \n",
					//  u, tau, exactrhs, abs(b[row] - exactrhs));
				}


				row++;
				if(row == rownum){
					break;
				}
			}

		}
	}

	

	if(row < rownum){
		fprintf(stderr, "Not enougth equations at (%i,%i,%i)\n", index[0],index[1],index[2]);
		SetMatrix(d2u,(double)NAN, 2, 2);
		SetMatrix(du, (double)NAN, 2);
	}

	vector<double> x = gesolve(A, b);
	u = x[0];

	for( int m = 0; m < 3; m ++){
		du[m] = x[DuCol(m)];
		for(int n = m; n < 3; n++){
				d2u[m][n] = x[D2uCol(m,n)];	
		}
	}

	if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
		printMat(10,10,A);
		printMat(10,b);
		PrintError(10, "solved", x.data(), "exact", xexact, "err");
	}
	
}


void CheckErrGeneral(double*** psi_real, double*** S, double*** psi_ghost, double*** sign_surf, PBData &pb, GridData &grid){
	cout<<"[CheckErrGeneral]"<<endl;
	int nancount = 0; // number of interface points fail
	
	double max_du_err = 0.0;
	Index max_du_idx;
	vector<double> du_err_all;
	

	double max_u_err = 0.0;
	Index max_u_idx;
	vector<double> u_err_all;


	for(int i = 0; i <= grid.nx[0]; i ++){
	for(int j = 0; j <= grid.nx[0]; j ++){
	for(int k = 0; k <= grid.nx[0]; k ++){
		Index index = {i, j, k};

		if(nearinterface(index, S, grid)){
			// For interface points
			
			// get exact u at grid point
			double exact_u = getu(index.data(), 0, 0, 0, GetSign(S, index), grid);

			double cu;
			double cdu[3];
			double** cd2u = matrix(2, 2);
			Vector3d cdu_itf; // computed du at interface


			if(IsGhost(index, sign_surf, S)){
				ComputeIcimLinSolve(cu, cdu, cd2u, psi_ghost, sign_surf, S, index, pb, grid);	
			}else{
				exit(1); // unfinished
			}

			// u error
			double u_err = abs(cu - exact_u);
			u_err_all.push_back(u_err);
			
			if(u_err > max_u_err){
				max_u_err = u_err;
				max_u_idx = index;
			}	
			

			for(int rstar = 0; rstar < 3; rstar ++){
				for(int sstar = -1; sstar <=1; sstar+=2){
					// get direction of interface
					Index rindex = UnitIncrement(index, rstar, sstar);
					if(!SameSide(S, index, rindex )){
						double alpha;
						double tangent[3], normal[3];
						Vector3d exact_Du = {0,0,0};
						
						getinterfaceinfo(alpha,tangent,normal,S,index.data(),rindex.data(),grid);

						// get excat Du at itnerface
						for(int r = 0; r < 3; r ++){
							exact_Du[r] = getDu(index.data(), r, rstar, sstar, alpha, GetSign(S, index), grid);
						}

						Vector3d h = {0.0};
						h[rstar] = alpha * grid.dx[rstar] * sstar;
						for(int i = 0; i < 3; i ++){
							cdu_itf[i] = TaylorExpand(cdu[i], cd2u[i], nullptr, h.data());
						}
					
						// du error
						for(int r = 0; r < 3; r ++){
							du_err_all.push_back(abs(cdu_itf[r] - exact_Du[r]));
						}

						double du_err = NormInf(Minus(cdu_itf, exact_Du));

						if(du_err > max_du_err){
							max_du_err = du_err;
							max_du_idx = index;
						}
					}
				}
			}

			free_matrix(cd2u, 2, 2);

		}

	}
	}
	}

	printf("u max err = %.8f at %i %i %i\n", max_u_err, max_u_idx[0], max_u_idx[1], max_u_idx[2]);
	printf("u rmse = %.8f \n", rmse(u_err_all));
	printf("Du max err = %.8f at %i %i %i\n", max_du_err, max_du_idx[0], max_du_idx[1], max_du_idx[2]);
	printf("Du rmse = %.8f \n", rmse(du_err_all));

}


void reconstruct(double*** psi_real, double*** S, double*** psi_ghost, double*** sign_surf, GridData& grid){
	for ( auto &it : D2uMap){
		Index ind;
		ind2sub(ind.data(), it.first, grid.nx, 3);
		RecoverD2u(ind, psi_ghost, grid);
	}


	for(int i = 0 ; i <= grid.nx[0]; i++){
		for(int j = 0 ; j <= grid.nx[1]; j++){
			for(int k = 0 ; k <= grid.nx[2]; k++){
				Index ind = {i,j,k};
				double uval = NAN;
				if (IsGhost(ind, S, sign_surf)){
					Vector3d dummy;
					uval = ComputeSimple(psi_ghost, sign_surf, S, ind, grid);
					if(isnan(uval)){
						ComputeByExtrapolate(uval, dummy, psi_ghost, sign_surf, S, ind, 0, 0, 0, grid);	
					}

				}else{
					// if ind is not not ghost state
					uval = psi_ghost[i][j][k];
				}

				if(isnan(uval)){
					cerr<<"fail to recontruct"<<endl;
					exit(1);
				}
				psi_real[i][j][k] = uval;
			}
		}
	}

}