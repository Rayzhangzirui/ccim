#include "catch.hpp"
#include "icim.h"
#include <limits>

extern int globcim;
extern double RADIUS;
extern int GRIDNUM;
extern int SURFOPT;
extern int argc;
extern char** argv;
extern double EPSILONP;
extern double EPSILONM;
extern double globGSsmooth;
extern int eindex[3];

int MAXITER = 10;
int DEPTH = 1;
double TOL = 1e-6;

void create_2sphere(double*** surf, GridData &grid, double r1, array<double,3> c1, double r2, array<double,3> c2){
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
				surf[i][j][k] = min(dist(x,c1) - r1, dist(x,c2) - r2);
			}
		}
	}
}

// center at c, axis along v, radius r
void create_cylinder(double*** surf, GridData &grid, double r, array<double,3> c, array<double,3> v){
	array<double,3> n = Div(v, Norm2(v));
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
				double dist = Norm2( Minus( Minus(c,x),  Mult(Dot(Minus(c,x), n) , n) ) );
				surf[i][j][k] = min( dist - r,  surf[i][j][k]);
			}
		}
	}
}



//check type function
TEST_CASE("type1"){
	double a = 1.0;
	GridData grid;
	init_grid(grid, 2 , a);
	
	double ***S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	double ***sign_surf = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	
	PBData pb;
	init_PBData(pb, grid, EPSILONP, EPSILONM);
	
	// fig5 (a)
	create_2sphere(S, grid, 1.5, array<double,3> {-1,1,1},1.5, array<double,3> {1,-1,-1});
	SignSurf(sign_surf, S, grid);
	print_surf(sign_surf, array<int,3> {1,1,1}, 1);

	// print Indicator
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				array<int,3> ind = {i,j,k};
				array<int,3> g = Indicator(ind, sign_surf, grid);
				print(ind);
				cout<<"-";
				print(g);
				int type = GetType(ind, sign_surf, grid);
				printf(" - %d",type);
				cout<<endl;

				// middle point shoudl be type 1
				if(ind[0]==1 && ind[1]==1 && ind[2]==1){
					CHECK(type==1);
				}else{
					CHECK(type==4);
				}
			}
		}
	}

	// flipping procedure
	flip(sign_surf, S, pb, grid, MAXITER, DEPTH);
	print_surf(sign_surf, array<int,3> {1,1,1}, 1);

}


TEST_CASE("type2a"){
	double a = 1.0;
	GridData grid;
	init_grid(grid, 2 , a);
	
	double*** S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	double*** sign_surf = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	
	PBData pb;
	init_PBData(pb, grid, EPSILONP, EPSILONM);

	// fig6 (a)
	cout<<"fig6 (a)"<<endl;
	SetField(S, grid, std::numeric_limits<double>::infinity());
	create_cylinder(S, grid, 0.5, array<double,3> {-1,-1,-1}, array<double,3>{0,0,1} );
	create_cylinder(S, grid, 1.1, array<double,3> {1,1,-1}, array<double,3>{0,0,1} );


	SignSurf(sign_surf, S, grid);
	print_surf(sign_surf, array<int,3> {1,1,1}, 1);

	// print Indicator
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				array<int,3> ind = {i,j,k};
				array<int,3> g = Indicator(ind, sign_surf, grid);
				print(ind);
				cout<<"-";
				print(g);
				int type = GetType(ind, sign_surf, grid);
				printf(" - %d", type);
				cout<<endl;

				// middle point shoudl be type 2
				if(ind[0]==1 && ind[1]==1 && ind[2]==1){
					CHECK(type==2);
				}else{
					CHECK(type==4);
				}
			}
		}
	}


	// flipping procedure
	
	flip(sign_surf, S, pb, grid, MAXITER, DEPTH);
	print_surf(sign_surf, array<int,3> {1,1,1}, 1);
}

TEST_CASE("type2f"){
	double a = 1.0;
	GridData grid;
	init_grid(grid, 2 , a);
	
	double*** S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	double*** sign_surf = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);

	PBData pb;
	init_PBData(pb, grid, EPSILONP, EPSILONM);

	// fig6 (f)
	cout<<"fig6 (f)"<<endl;
	SetField(S, grid, std::numeric_limits<double>::infinity());
	create_cylinder(S, grid, 0.5, array<double,3> {-1,-1,-1}, array<double,3>{0,0,1} );
	create_cylinder(S, grid, 0.51, array<double,3> {1,-0.5,-1}, array<double,3>{0,0,1} );
	create_cylinder(S, grid, 0.51, array<double,3> {-0.5,1,-1}, array<double,3>{0,0,1} );


	SignSurf(sign_surf, S, grid);
	print_surf(sign_surf, array<int,3> {1,1,1}, 1);

	// print Indicator
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				array<int,3> ind = {i,j,k};
				array<int,3> g = Indicator(ind, sign_surf, grid);
				print(ind);
				cout<<" - ";
				print(g);
				int type = GetType(ind, sign_surf, grid);
				printf(" - %d", type);
				cout<<endl;

				// middle point shoudl be type 2
				if(ind[0]==1 && ind[1]==1 && ind[2]==1){
					CHECK(type==2);
				}else{
					CHECK(type==4);
				}
			}
		}
	}


	// flipping procedure
	
	
	flip(sign_surf, S, pb, grid, MAXITER, DEPTH);
	print_surf(sign_surf, array<int,3> {1,1,1}, 1);

}

TEST_CASE("flip"){

	cmdLine(argc, argv);// can not take in arguments, only print

	double tol = 1.0e-14; // tol for pertubation. tol for solver in grid
	double gridbound = 1.0;

	int count_flip[4] = {0};
	int count_noflip[4] = {0};
	int no_fix = 0;

	GridData grid;
	init_grid(grid, GRIDNUM, gridbound);
	
	PBData pb;
	init_PBData(pb, grid, EPSILONP, EPSILONM);

	double ***S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	init_surf_perturb(S, grid.radius, grid, SURFOPT, pb, tol);

	double*** sign_surf = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	SignSurf(sign_surf, S, grid);
	
	
	flip(sign_surf, S, pb, grid, MAXITER, DEPTH);

	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				array<int,3> ind = {i,j,k};
				int type = GetType(ind, sign_surf, grid);

				int status_flip = getstatus3(sign_surf, ind.data(), grid);

				int status_noflip = getstatus3(S, ind.data(), grid);
				
				if(status_noflip == 0){ 
				//boundary
					count_noflip[0]++;
				}

				if(status_noflip == 1){
				//interior
					count_noflip[1]++;
				}

				if(status_noflip == 2){
					// usual cim2
					CHECK(IsUsualCim2(ind, S, grid)==true);
					count_noflip[2]++;
				}

				if(status_noflip == 3){
					// exceptional
					count_noflip[3]++;
				}



				if(status_flip == 0){ 
				//boundary
					CHECK(type==4);
					count_flip[0]++;
				}

				if(status_flip == 1){
				//interior
					CHECK(type == 0);
					count_flip[1]++;
				}

				if(status_flip == 2){ // after flipping, should be cim2 points
					//cim2 after flip
					CHECK(IsIcim2(ind, sign_surf, grid)==true);
					

					count_flip[2]++;
				}

				if(status_flip == 3){
					// cim2 after flip still exceptional

					if (IsIcim2(ind, S, grid) == false){
						printf("not icim2 at (%i,%i,%i) type = %i\n", ind[0],ind[1],ind[2], type);
						print_surf(sign_surf,ind,2);
						cout<<endl;
						no_fix++;
					}else{
						count_flip[3]++;
					}

					
				}
			}
		}
	}


	cout<<"\nbefore flip"<<endl;
	cout<<"boundary points = " << count_noflip[0] <<endl;
	cout<<"interior points = " << count_noflip[1] <<endl;
	cout<<"cim2 points = " << count_noflip[2] <<endl;
	cout<<"exceptional points = " << count_noflip[3] <<endl;

	cout<<"\nafter flip"<<endl;
	cout<<"boundary points = " << count_flip[0] <<endl;
	cout<<"interior points = " << count_flip[1] <<endl;
	cout<<"usual cim2 points = " << count_flip[2] <<endl;
	cout<<"icim2 points = " << count_flip[3] <<endl;
	cout<<"no fix = " << no_fix <<endl;

	CHECK(no_fix == 0); // should have no exceptional after flipping
	

	free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
	free_matrix(sign_surf,grid.nx[0],grid.nx[1],grid.nx[2]);


}

// test local truncation error
TEST_CASE("get"){

	cmdLine(argc, argv);// can not take in arguments, only print

	double tol = 1.0e-12; // tol for pertubation. tol for solver in grid
	double gridbound = 1.0;

	GridData grid;
	init_grid(grid, GRIDNUM, gridbound);
	
	PBData pb;
	init_PBData(pb, grid, EPSILONP, EPSILONM);

	double ***S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	init_surf_perturb(S, grid.radius, grid, SURFOPT, pb, tol);


	// approximate D2u
	int mid = 2;
	int N = 2 * mid;

	double***** D2ucoef = new double ****[grid.dim];
	double*** D2uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
	double tempuxcoef[3] = {0,0,0};
	for (int r = 0; r < grid.dim; r++)
	{
		D2ucoef[r] = new double ***[grid.dim];
		for (int s = 0; s < grid.dim; s++)
			D2ucoef[r][s] = matrix(N,N,N);
	}

	double err;
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				Index ind = {i,j,k};
				// work on interface pointts
				if(nearinterface(ind,  S, grid)){

					// for each plane, check D2u
					for (int k = 0; k < 3; k++){
						for(int l = k + 1; l < 3; l++){
							
							double exactD2u = getD2u(ind.data(),k,l,0,0,0,GetSign(S,ind),grid);
							// usual stencil
							bool has_usual = HasUsualMixD2u(D2ucoef[k][l], ind, k, l, S,  grid);
							double apprxD2u_usual = evalcoef( 0, D2ucoef[k][l], tempuxcoef, D2uxxcoef[k][l], ind.data(), 0, 0, 0, mid, S, grid);
							err = exactD2u-apprxD2u_usual;

							if(has_usual){
								if(abs(err)>TOL){
									fprintf(stderr, "D2u[%i][%i] at (%i,%i,%i)"
										" exactD2u = %f , apprxD2u_usual = %f, err = %f\n", k, l, ind[0],ind[1],ind[2],exactD2u,apprxD2u_usual,err);
									exit(1);
								}
								
							}

							// couple with second derivitive
							bool has_couple = HasCoupleMixD2u(D2ucoef[k][l], D2uxxcoef[k][l], ind, k, l, S, grid);
							double apprxD2u_couple = evalcoef( 0, D2ucoef[k][l], tempuxcoef, D2uxxcoef[k][l], ind.data(), 0, 0, 0, mid, S, grid);
							err = exactD2u - apprxD2u_couple;
							
							if(has_couple){
								if(abs(err)>TOL){
									fprintf(stderr, "D2u[%i][%i] at (%i,%i,%i)"
										" exactD2u = %f , apprxD2u_couple = %f, err = %f\n", k, l, ind[0],ind[1],ind[2],exactD2u,apprxD2u_couple,err);
									exit(1);
								}
								
							}

							// shift to nbr
							bool has_shift = HasShiftMixD2u( D2ucoef[k][l], ind, k, l, S, grid);
							
							double apprxD2u_shift = evalcoef( 0, D2ucoef[k][l], tempuxcoef, D2uxxcoef[k][l], ind.data(), 0, 0, 0, mid, S, grid);
							err = exactD2u - apprxD2u_shift;
							
							if(has_couple){
								if(abs(err)>TOL){
									fprintf(stderr, "D2u[%i][%i] at (%i,%i,%i)"
										" exactD2u = %f , apprxD2u_shift = %f, err = %f\n", k, l, ind[0],ind[1],ind[2],exactD2u,apprxD2u_shift,err);
									exit(1);
								}
								
							}

							if( (!has_couple) && (!has_usual) && (!has_shift) ){
								fprintf(stderr, "no D2u[%i][%i] at (%i,%i,%i)\n", k, l, ind[0],ind[1],ind[2]);
								printPlane(S, grid, ind.data(), k, l, 2);
							}

							
						}
					}

					double**** D1ucoef = new double ***[grid.dim];
					for (int r = 0; r < grid.dim; r++)
						D1ucoef[r] = matrix(N,N,N);
					double** D1uxxcoef = matrix(grid.dim-1,grid.dim-1);
					double dummy_D1uxcoef[3] = {0};
					double exact_Du[3];

					for(int m = 0; m < 3; m++){						
						exact_Du[m] = getDu(ind.data(),m,0,0,0,GetSign(S,ind),grid);
					}
					// check Du at gridpoint
					GetIcimDuGridpoint(D1ucoef, D1uxxcoef, S, ind.data(), grid);
				
					// get exact_Du
					for(int m = 0; m < 3; m++){						
						double approx = evalcoef(0.0, D1ucoef[m],dummy_D1uxcoef,D1uxxcoef[m],ind.data(),0,0,0.0,mid,S,grid);
						if(abs(err)>TOL){
							printf("(%i, %i, %i), Du[%i] = %f, exact = %f, err =  %f\n", ind[0], ind[1], ind[2], m, approx, exact_Du, err);
						}
						
					}
					

					// for each direction check D1u at interface
					for (int r = 0; r < 3; r ++){
						for (int s = -1; s <= 2; s += 2){
							Index rindex = UnitIncrement(ind, r, s);
							if (!SameSide(S, ind, rindex)){
								
								double alpha;
								double tangent[3], normal[3];
								getinterfaceinfo(alpha,tangent,normal,S,ind.data(),rindex.data(),grid);

								// Get coeff
								GetIcimDu(D1ucoef, D1uxxcoef, S, ind.data(), r, s, alpha, grid);

								// check Du each component finite difference 
								for(int m = 0; m < 3; m++){
									double thesign = GetSign(S,ind);
									double approx = evalcoef(0.0, D1ucoef[m],dummy_D1uxcoef,D1uxxcoef[m],ind.data(),0,0,0.0,mid,S,grid);
									double exact_Du = getDu(ind.data(),m,r,s,alpha,thesign,grid);
									
									double err_apprx = exact_Du - approx;

									if(abs(err_apprx)>TOL){
										printf("(%i, %i, %i), rstar = %i, sstar = %i: Du[%i] = %f, exact = %f, err_apprx =  %f\n",
										 ind[0], ind[1], ind[2], r, s, m, approx, exact_Du, err_apprx);
										cout<<"u coef"<<endl;
										PrintCoef(D1ucoef[m],4);
										cout<<"uxx coef"<<endl;
										PrintCoef(D1uxxcoef[m],2);
										exit(1);
									}
								}
							}
						}
					}

					
				}
			}
		}
	}


	
	free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
	
}


// check reconstruction
TEST_CASE("compute"){

	cmdLine(argc, argv);// can not take in arguments, only print

	eindex[0] = 7;
	eindex[1] = 9;
	eindex[2] = 9;

	double tol = 1.0e-14; // tol for pertubation. tol for solver in grid
	double gridbound = 1.0;


	GridData grid;
	init_grid(grid, GRIDNUM, gridbound);
	
	PBData pb;
	init_PBData(pb, grid, EPSILONP, EPSILONM);

	double ***S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	init_surf_perturb(S, grid.radius, grid, SURFOPT, pb, tol);
	
	double*** u_exact = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				int index[3] = {i, j, k};
				u_exact[i][j][k] = getu(index, 0, 0, 0, GetSign(S,index),grid);
			}		
		}
	}


	double err;
	double** D2u = matrix(2,2);
	double Du[3];
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				Index ind = {i,j,k};
				// work on interface pointts
				if(!nearinterface(ind,  S, grid)){
					continue;
				}

				// mix D2u at grid point
				ComputeMixD2u(D2u, u_exact, u_exact, S, S, ind, grid);
				// for each plane, check D2u
				for (int m = 0; m < 3; m++){
					for(int n = m + 1; n < 3; n++){

						double exact_uxx = getD2u(ind.data(),m,n,0,0,0,GetSign(S,ind),grid);
						double err = D2u[m][n] - exact_uxx;

						if(abs(err)>TOL){
							printf("(%i, %i, %i) ",ind[0],ind[1],ind[2]);
							printf("D2u[%i][%i] compute = %f, exact = %f, err = %f\n", m, n, D2u[m][n], exact_uxx, abs(D2u[m][n] - exact_uxx) );	
							exit(1);
						}

					}
				}

				// du at gridpoint
				ComputeDuGridpoint(Du, u_exact, u_exact, S, S, ind, grid);

				double exact_du_gp[3];
				for(int r = 0; r < 3; r ++){
					exact_du_gp[r] = getDu(ind.data(), r, 0, 0, 0, GetSign(S, ind), grid);
					double err = exact_du_gp[r] - Du[r];
					
					if(abs(err)>TOL){
						printf("(%i, %i, %i) ",ind[0],ind[1],ind[2]);
						printf("Du[%i] compute = %f, exact = %f, err = %f\n", r, Du[r], exact_du_gp[r], abs(err) );	
						exit(1);
					}
					
				}


				// Du and u at interface
				for (int r = 0; r < 3; r ++){
					for (int s = -1; s <= 2; s += 2){
						Index rindex = UnitIncrement(ind, r, s);
						if (SameSide(S, ind, rindex)){
							continue;
						}
							
						double alpha;
						double tangent[3], normal[3];
						getinterfaceinfo(alpha,tangent,normal,S,ind.data(),rindex.data(),grid);

						Vector3d duitf = ComputeIcimDu(u_exact, u_exact, S, S, ind, r, s, alpha, grid);

						double exact_du_itf[3];
						// check each component
						for(int j = 0; j < 3; j++){
							exact_du_itf[j] = getDu(ind.data(),j,r,s,alpha,GetSign(S,ind),grid);
							double err = duitf[j] - exact_du_itf[j];
							if(abs(err)>TOL){
								printf("(%i, %i, %i) ",ind[0],ind[1],ind[2]);
								printf("Du[%i](xhat) exact = %f, compute = %f,  err = %f.\n", j, exact_du_itf[j], Du[j], err );	
								exit(1);
							}
						}


						double uitf = ComputeIcimU(u_exact, u_exact, S, S, ind, r, s, alpha, grid);


						// check each component
						
						double exact_u = getu(ind.data(),r,s,alpha,GetSign(S,ind),grid);
						double err = uitf - exact_u;
						if(abs(err)>TOL){
							printf("(%i, %i, %i) ",ind[0],ind[1],ind[2]);
							printf("u(xhat) exact = %f, compute = %f,  err = %f.\n", exact_u, uitf, err );	
							exit(1);
						}
						


					}
				}
				// end of check du and u at interface

			}
		}
	}


	
	free_matrix(S,grid.nx[0],grid.nx[1],grid.nx[2]);
	
}


TEST_CASE("icim"){
	eindex[0] = 10;
	eindex[1] = 12;
	eindex[2] = 10;
	globcim = 2;
	cmdLine(argc, argv);
	

	double tol = 1.0e-14; // tol for pertubation. tol for solver in grid
	double gridbound = 1.0;

	GridData grid;
	init_grid(grid, GRIDNUM, gridbound);
	
	PBData pb;
	init_PBData(pb, grid, EPSILONP, EPSILONM);

	double ***S = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	init_surf_perturb(S, grid.radius, grid, SURFOPT, pb, tol);

	double*** sign_surf = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	SignSurf(sign_surf, S, grid);

	
	flip(sign_surf, S, pb, grid, MAXITER, DEPTH , true);

	TempStruct tmp;//sparse linear solver 
	init_TempStruct(tmp, grid);

	StorageStruct *Dusmall;
	int smallsize = 0;
	// linearsystem5(tmp.A, tmp.b, Dusmall, smallsize, S, pb, grid);
	linearsystem_icim(tmp.A, tmp.b, S, sign_surf, pb, grid);


	// solve linear system with prerightBICGSTAB
	ILUsmall(tmp.M,tmp.A,grid,sign_surf,pb);
	prerightBICGSTABsmall(pb.psi,tmp.A,tmp.b,tmp.M,grid, grid.nx[0]*grid.nx[1]*grid.nx[2],1.0e-14, sign_surf, pb, tmp);
	gaussseidelsmall(pb.psi,tmp.A,tmp.b,grid,globGSsmooth, sign_surf, pb);

	// check result
	cout<<"Ghost state u error"<<endl;
	checkanswer(pb.psi, sign_surf, pb, grid); // check with flipped surface

	double *** psi_true = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
	reconstruct(psi_true, S, pb.psi,sign_surf, grid);

	cout<<"Actual u error"<<endl;
	checkanswer(psi_true, S, pb, grid); // check original surface	

	CheckIcimDu(pb.psi, psi_true, sign_surf, S, pb, grid);
	
}