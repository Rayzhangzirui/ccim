// #define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
// #define CATCH_CONFIG_RUNNER
// #define CATCH_CONFIG_ENABLE_BENCHMARKING
#include "catch.hpp"
#include "helper.h"
#include "icim.h"
#include "tryvn.h"
#include <random>

using namespace std;

extern int globtestnum;
extern string PROTEINFILE;
extern string PROJECTDIR;
extern int GRIDNUM;
// write surface
TEST_CASE("write_surf"){
	GridData grid;
	int GRIDNUM = 50;
	double a = 1.0;
	init_grid(grid, GRIDNUM , a);
	double radius = 0.5;
	// create initial surface
	double ***S;
	S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	for(int opt = 11; opt<=16; opt++){
		init_surf(S, radius, grid, opt);
		write_field("initialsurface"+to_string(opt)+".txt", S, grid);	
	}	
	int opt = 1;
	init_surf(S, radius, grid, opt);
	write_field("initialsurface"+to_string(opt)+".txt", S, grid);	
}

// creat protein surface with smoothing method in shu's paper
TEST_CASE("papersurf","[psurf]"){
	GridData grid;
	GRIDNUM = 100;
	double a = 1.0;
	init_grid(grid, GRIDNUM , a);
	
	double ***S;
	S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	init_surf_protein_paper(S, grid);
	string fname("surf_protein_paper");
	string filepath = PROJECTDIR + fname + "_gd" + to_string(GRIDNUM) + "_" + PROTEINFILE;
	write_field(filepath, S, grid);
	free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
}

// creat protein surface with heatsmooth
TEST_CASE("heatsurf","[psurf]"){
	GridData grid;
	GRIDNUM = 100;
	double a = 1.0;
	init_grid(grid, GRIDNUM , a);
		
	for(int heatstep = 0; heatstep <=2; heatstep+=1){
		double ***S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
		init_surf_protein_dist(S, grid, heatstep);
		string fname("surf_protein_dist");
		string filepath = PROJECTDIR + fname + "_gd" + to_string(GRIDNUM) + "_" + PROTEINFILE;
		write_field(filepath, S, grid);	
		free_matrix(S, grid.nx[0],grid.nx[1],grid.nx[2]);
	}
}

// investigate how cim1 points occur in protein surface
TEST_CASE("cim1"){
	GridData grid;
	int GRIDNUM = 185;
	double a = 1.0;
	init_grid(grid, GRIDNUM , a);
	// create initial surface
	double ***S;
	S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	init_surf_protein_paper(S, grid);
	int index[3] = {70,79,137};
	int status = getstatus5debug(S,index,grid);
	cout<<"status = " <<status<<endl;
	print_surf(S,index,2);
}


// investigate how cim1 points occur in protein surface
TEST_CASE("banana"){
	GridData grid;
	int GRIDNUM = 110;
	double a = 1.0;
	init_grid(grid, GRIDNUM , a);
	// create initial surface
	double ***S;
	S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	
	double radius = 0.5;
	int opt = 13;
	init_surf(S, radius, grid, opt);
	int index[3] = {18,54,11};
	int status = getstatus5debug(S,index,grid);
	cout<<"status = " <<status<<endl;
	print_surf(S,index,2);
}

// test condition number estimation with small matrix
TEST_CASE("small matrix"){
	int dim = 3;
	vector<vector<double>> A(dim, vector<double> (dim, 0));
	
	for(int i = 0; i < dim; i++){
		A[i][i] = 1.0;
	}
	CHECK(Approx(1.0)==condest(A));
	
	for(int i = 0; i < dim; i++){
		A[i][i] = 1.0 + i;
	}
	CHECK(Approx(1.0*dim) == condest(A));

	for(int i = 0; i < dim; i++){
		A[i][i] = 0.0;
	}
	A[0][0] = 1.0;A[0][1] = 2.0;A[0][2] = 3.0;
				  A[1][1] = 4.0;A[1][2] = 5.0;
								A[2][2] = 6.0;
	CHECK(Approx(14.0)==condest(A));
}

TEST_CASE("D2u"){
	// in yessk2
	cout << "yessk2"<<endl;
	int tempsk2[4];
	bool foundone = false;
	for (int m = 2; m >= 1 && !foundone; m--)
		for (int n = 2; n >= 1 && !foundone; n--)
			for (tempsk2[0] = -1,tempsk2[1] = tempsk2[0]+m; tempsk2[1] <= 1; (tempsk2[0])++,tempsk2[1] = tempsk2[0]+m)
				for (tempsk2[2] = -1,tempsk2[3] = tempsk2[2]+n; tempsk2[3] <= 1; (tempsk2[2])++,tempsk2[3] = tempsk2[2]+n)
					cout<<tempsk2[0]<<" "<<tempsk2[1]<<" "<<tempsk2[2]<<" "<<tempsk2[3]<<endl;
	

	cout << "yescim5D2"<<endl;
	int ntheta = 8;
	double dtheta = 2.0*M_PI/ntheta;
	double theta0 = -3.0*M_PI/4.0;// start from bottom left corner
	int index[3] = {0,0,0};
	int tindex[3] = {0,0,0};;
	int offset[2][2];
	int m = 0, n = 1;
	for (int r = 0; r < ntheta; r++)
	{
		double theta = theta0+r*dtheta;
		double themax = max(fabs(cos(theta)),fabs(sin(theta)));
		offset[0][0] = round(cos(theta)/themax);
		offset[0][1] = round(sin(theta)/themax);
		tindex[m] = index[m]+offset[0][0];
		tindex[n] = index[n]+offset[0][1];
		printf("tindex (%d,%d,%d)\n",tindex[0],tindex[1],tindex[2]);
	}




}

// test condition number estimation with random
TEST_CASE("random matrix"){
	int tn = 10000;
	int dim = 6;
	ofstream outfile("test_mx.txt");
	for(int i = 0; i<tn; i++){
		vector<vector<double>> A = randmx(dim,dim);
		vector<double> b = randmx(dim);
		vector<double> x = gesolve(A,b);
		//write output
		for(int j = 0; j < dim*dim; j++){
			outfile<<std::setprecision(15) <<A[j/dim][j%dim]<<" ";
		}
		for(int j = 0; j< dim; j++){
			outfile<<b[j]<<" ";
		}
		for(int j = 0; j< dim; j++){
			outfile<<x[j]<<" ";
		}
		outfile<<norm1(A)<<" "<<condest(A)<<endl;
	}
	
}





// investigate different approximation of D2u
TEST_CASE("jump"){
	
	
	// banana shape
	int GRIDNUM = 110;
	int opt = 13;
	int index[3] = {18,54,11};
	
	// peanut shape
	// int GRIDNUM = 70;
	// int opt = 15;
	// int index[3] = {36,36,17};



	GridData grid;
	globtestnum = 11;
	
	init_grid(grid, GRIDNUM , 1.0);
	double ***S;
	S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	
	
	init_surf(S, 0.0, grid, opt);
	


	
	double gamma[3][2],alpha[3][2],tangent[3], normal[3];

	double diffscheme[5][4] =  {{-1, 1, -1, 1},
								{-1, 0, -1, 0},
								{-1, 0,  0, 1},
								{ 0, 1, -1, 0},
								{ 0, 1,  0, 1}};

	int rindex[3];
	copy(index, index+3, rindex);

	int mid = 2, N = 2*mid;

	double **D2u, *****D2ucoef, ***D2uxcoef, ***D2uxxcoef, **D2jumpuxxcoef;
	D2u = matrix(grid.dim-1,grid.dim-1);
	D2ucoef = new double ****[grid.dim];
	for (int r = 0; r < grid.dim; r++){
		D2ucoef[r] = new double ***[grid.dim];
		for (int s = 0; s < grid.dim; s++)
			D2ucoef[r][s] = matrix(N,N,N);
	}
	D2uxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1); //3x3x3
	D2uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);//3x3x3
	D2jumpuxxcoef = matrix(grid.dim-1,grid.dim-1);//3x3 

	// calculate gamma
	for (int r = 0; r < grid.dim; r++){
		for (int s = -1; s <= 1; s += 2){
			rindex[r] = min(max(index[r]+s,0),grid.nx[r]);//if rindex (nbr of index) goes out of boundary, set as index
			if (evalarray(S,index)*evalarray(S,rindex) < 0.0)
				gamma[r][(s+1)/2] = 1;
			else
				gamma[r][(s+1)/2] = 0;
		}
		rindex[r] = index[r];
	}

	for (int m = 0; m < grid.dim; m++){
		for (int n = m+1; n < grid.dim; n++){
			printf("Plane %d %d surface\n",m,n);
			printPlane(S,grid,index,m,n,2);
		}
	}

	double thesign = (evalarray(S,index) < 0.0)? -1.0: 1.0;

	// finite difference to approximate cross dierivative
	for (int r = 0; r < grid.dim; r++){
		for (int sk = -1; sk <= 1; sk += 2){
			copy(index, index+3, rindex);
			if (gamma[r][(sk+1)/2] == 1){
				
				rindex[r] = index[r]+sk;

				getinterfaceinfo(alpha[r][(sk+1)/2],tangent,normal,S,index,rindex,grid);
				double beta = 1.0-alpha[r][(sk+1)/2];

				double x[3];
				sub2coord(x,index,grid);
				x[r] += sk*alpha[r][(sk+1)/2]*grid.dx[r];
				printf("Interface in r = %d， sk = %d at (%f,%f,%f), alpha = %f\n", r,sk,x[0],x[1],x[2],alpha[r][(sk+1)/2]);

				//for each plane D2u[m][n]
				for (int m = 0; m < grid.dim; m++){
					for (int n = m+1; n < grid.dim; n++){
						

						double D2ueItf = getD2u(index,m,n,r,sk,alpha[r][(sk+1)/2],thesign,grid);//exact D2u at interface
						double D2ueGp = getD2u(index,m,n,r,sk,0,thesign,grid);//exact D2u at grid point

						printf("	D2u[%d][%d] exact at grid point = %f, exact at interface = %f\n",m,n,D2ueGp,D2ueItf);
						
						// for each difference scheme
						for (int k = 0; k < 5; k ++){
							int sk2[4];
							
							copy(diffscheme[k],diffscheme[k]+4,sk2);

							int nindex[3] = {N,N,N};
							int sindex[3] = {mid,mid,mid};
							getD2(D2ucoef[m][n],m,n,sk2,sindex,nindex,grid);// approximate D2u in terms of u-value only

							for (int t = 0; t < grid.dim; t++){
								D2uxcoef[m][n][t] = 0.0;
								D2uxxcoef[m][n][t] = 0.0;
							}
							D2jumpuxxcoef[m][n] = 0.0;
							D2u[m][n] = 0;

							double D2uApx = evalcoef(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],index,0,0,0.0,mid,thesign,grid);//Apprx D2u at grid point
							
							printf("		cim3 sk2 = (%d,%d,%d,%d), D2u apprx = %f, error = %f\n",sk2[0], sk2[1], sk2[2], sk2[3], D2uApx, D2uApx - D2ueItf);
						}
						
						//use cim5
						if(getcim5D2(D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],m,n,index,mid,S,grid)){
							D2u[m][n] = 0;
							D2jumpuxxcoef[m][n] = 0.0;
							
							double D2uApx = evalcoef(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],index,0,0,0.0,mid,thesign,grid);//Apprx D2u at grid point
							printf("		cim5 D2u apprx = %f, error = %f\n", D2uApx, D2uApx - D2ueItf);
						}else{
							printf("		No cim5 scheme\n");
						}

						// get all cim5 scheme
						vector<double***> D2ucoefvec;
						vector<double*> D2uxcoefvec, D2uxxcoefvec;
						vector<vector<int>> offsetvec;
						bool yes = yescim5D2All(D2ucoefvec,D2uxcoefvec,D2uxxcoefvec,offsetvec,m,n,index,mid,S,grid);
						D2u[m][n] = 0;
						if(yes){
							for(int i = 0; i < D2ucoefvec.size(); i++){
								double D2uApx = evalcoef(D2u[m][n], D2ucoefvec[i],D2uxcoefvec[i],D2uxxcoefvec[i],index,0,0,0.0,mid,thesign,grid);
								printf("		cim5 D2u sk2 = (%d,%d,%d,%d) apprx = %f, error = %f\n",offsetvec[i][0],offsetvec[i][1],offsetvec[i][2],offsetvec[i][3], D2uApx, D2uApx - D2ueItf);
							}
						}


						// use the other side
						printPlane(S,grid,rindex,m,n,2);
						int sk2[4];
						if (getsk2(sk2,m,n,rindex,S,grid)){						
							int nindex[3] = {N,N,N};
							int sindex[3] = {mid,mid,mid};

							sindex[r] = mid+sk;
							getD2(D2ucoef[m][n],m,n,sk2,sindex,nindex,grid);
							
							

							for (int t = 0; t < grid.dim; t++){
								D2uxcoef[m][n][t] = 0.0;
								D2uxxcoef[m][n][t] = 0.0;
							}
							D2jumpuxxcoef[m][n] = 0.0;
							D2u[m][n] = 0;
							D2jumpuxxcoef[m][n] = thesign;						
							
							double D2uApx = evalcoef(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],index,0,0,0.0,mid,thesign,grid);//Apprx D2u at grid point
							printf("		cim4 sk2 = (%d,%d,%d,%d), D2u apprx = %f, error = %f\n",sk2[0], sk2[1], sk2[2], sk2[3], D2uApx, D2uApx - D2ueItf);
						}else{
							printf("		No other side scheme\n");
						}
						
					}
				}

						
			}
		}
	}


	free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
}


