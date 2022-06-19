//max difference
#include "helper.h"
#include <string>
#include <sstream>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <random>

using namespace std;


int SURFOPT = 0;
double RADIUS = 0.5;
int runStep = 2000;
double globtime = 0.1;
double ETA = 30.0;
string PROTEINDIR = "/Users/zzirui/pbcim/protein_surf/";
string PROJECTDIR = "/Users/zzirui/pbcim/";
string PROTEINFILE = "mdm2.txt";

// defined in tryvn.C
extern char globbiasstat;
extern int globtesta;
extern int globsmall;
extern int GRIDNUM;
extern int globtestnum;
extern int globperturb;
extern int globcim;
extern char globallcim345;
extern int globlinsolve;
extern double tollinsolve;
extern char globcheck;
extern int globGSsmooth;
extern int globheatsmooth;
extern char globdtdx2;
extern double globheatsmoothtime;
extern int eindex[3];
extern double EPSILONP;
extern double EPSILONM;
extern char globdist;
extern bool globwritemx;
extern ofstream outfilemx;
extern char globintorder;
bool globgridperturb = false;

//initialize grid
void init_grid(GridData &grid, int nx, double a){
	

	grid.nx[0] = nx; //number of cell
	grid.nx[1] = grid.nx[0];
	grid.nx[2] = grid.nx[0];
	grid.a[0] = -a;
	grid.a[1] = -a;
	grid.a[2] = -a;
	grid.dx[0] = 2.0*fabs(grid.a[0])/grid.nx[0];
	grid.dx[1] = 2.0*fabs(grid.a[1])/grid.nx[1];
	grid.dx[2] = 2.0*fabs(grid.a[2])/grid.nx[2];
  grid.mindx = min(min(grid.dx[0],grid.dx[1]),grid.dx[2]);
	grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim);
	grid.tol = 1.0e-14;
	grid.t = 0.0;

  grid.radius = RADIUS; // current radius for sphere
  grid.radius0 = RADIUS; // initial radius for sphere


  if(globgridperturb){
    std::default_random_engine generator(std::random_device{}());
    std::uniform_real_distribution<double> distribution(0.0, grid.mindx);
    double rand_shift[3] = {distribution(generator), distribution(generator), distribution(generator)};
    grid.a[0] += rand_shift[0];
    grid.a[1] += rand_shift[1];
    grid.a[2] += rand_shift[2];
    printf("randomly perturb grid by (%f,%f,%f)\n",rand_shift[0],rand_shift[1],rand_shift[2]);
  }
}


// Initialize sparse matrix solver
void init_TempStruct(TempStruct &tmp, GridData &grid){
	tmp.Nfd = 10;
   tmp.fourd = new double ***[tmp.Nfd];
   tmp.fdstatus = new char[tmp.Nfd];
   for (int i = 0; i < tmp.Nfd; i++)
   {
      tmp.fourd[i] = NULL;
      tmp.fdstatus[i] = 0;
   }
   tmp.A = sparseelt2ptrmatrix(grid.nx[0],grid.nx[1],grid.nx[2]); // A is 3d array,each is a pointer to sparseelt
   if (globlinsolve == 1){
      tmp.M = sparseelt2ptrmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   }
   tmp.b = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
}

void init_PBData(PBData &pb, GridData &grid, double epsp, double epsm){
	pb.epsilonm = epsm;
   pb.epsilonp = epsp;
   pb.gamma0 = 0.0;
   pb.psi = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   pb.dim = grid.dim;
   pb.N = 2;//no solute
   pb.x = matrix(pb.N-1,grid.dim-1);
   pb.x[0][0] = 0.1;
   pb.x[0][1] = 0.0;
   pb.x[0][2] = 0.0;
   pb.x[1][0] = -0.1;
   pb.x[1][1] = 0.0;
   pb.x[1][2] = 0.0;
   pb.Q = new double[pb.N];
   pb.c = new double[pb.N];
   pb.epsilone = new double[pb.N];
   pb.epsilonlj = new double[pb.N];
   pb.sigmalj = new double[pb.N];
   for (int i = 0; i < pb.N; i++)
   {
      pb.Q[i] = 0.0;//1.0;
      pb.c[i] = 0.0;//1.0;
      pb.epsilone[i] = 0.0;//1.0;
      pb.epsilonlj[i] = 0.0;//0.0159;
      pb.sigmalj[i] = 0.0;//3.653;
   }
   pb.beta = 1.0;
   pb.rho0 = 0.0333;
   pb.LJ = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);


   for(int i = 0; i <= grid.nx[0]; i++){
      for(int j = 0; j <= grid.nx[1]; j++){
         for(int k = 0; k <= grid.nx[2]; k++){
            pb.psi[i][j][k] = 0;
            pb.LJ[i][j][k] = 0;
         }
      }
   }

}



// create sphere surface
void create_sphere(double*** surf, GridData &grid, double radius, std::array<double,3> center){
	for(int i = 0; i <= grid.nx[0]; i++){
		for(int j = 0; j <= grid.nx[1]; j++){
			for(int k = 0; k <= grid.nx[2]; k++){
				array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
				surf[i][j][k] = dist(x,center) - radius;
			}
		}
	}
}

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


void init_surf(double ***S, double radius, GridData &grid, int opt){
   if(opt == 0){
      create_sphere(S,grid,radius,array<double,3>{0,0,0});   
   }else if (opt == 1){
      // strange torus from getinit
      double radii[3];
      radii[0] = 0.5;
      radii[1] = 0.25;
      radii[2] = 0.15;
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
               double value = 0;
               for (int r = 1; r < grid.dim; r++)
                  value += x[r]*x[r];
               value = sqrt(value)-radii[0];
               double theta = atan2(x[2],x[1]);
   // main one used in tests
               value = sqrt(value*value+x[0]*x[0])-(0.5*(radii[1]-radii[2])*sin(3.0*theta)+
                                              0.5*(radii[1]+radii[2]));
   //      value = sqrt(value*value+x[0]*x[0])-(0.5*(radii[1]-radii[2])*sin(theta)+
   //                                           0.5*(radii[1]+radii[2]));
   //      value = sqrt(value*value+x[0]*x[0])-radii[1];
               S[i][j][k] = value;
            }
         }
      }

   }
   else if (opt == 11){
      //ellipsoid
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
               S[i][j][k]= 2 * pow(x[0],2) + 3 * pow(x[1],2) + 6 * pow(x[2],2) - 1.3 * 1.3;
            }
         }
      }
   }else if (opt == 12){
      //donut from paper
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
               S[i][j][k]= -0.3 * 0.3 + pow(-0.6 + sqrt(pow(x[0],2) + pow(x[1],2)),2) + pow(x[2],2);
            }
         }
      }
   }else if (opt == 13){
      //banana
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
               S[i][j][k] = 1521 - 94*pow(6 + 7*x[0],2) + pow(6 + 7*x[0],4) + 3822*pow(x[1],2) + 2401*pow(x[1],4) - 4606*pow(x[2],2) + 4802*pow(x[1],2)*pow(x[2],2) + 3601.5*pow(x[2],4) + 98*pow(6 + 7*x[0],2)*(pow(x[1],2) + pow(x[2],2));
            }
         }
      }
   }else if (opt == 14){
      //8 balls
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
               double d = std::numeric_limits<double>::infinity();
               for(int n = 0; n <=7; n++){
                  array<double,3> c = {pow(-1,floor(n/4)) * 0.5, pow(-1,floor(n/2)) * 0.5, pow(-1,n) * 0.5};
                  d = min(d,dist(x,c)-0.3);
                  S[i][j][k] = d;
               }
            }
         }
      }
   }else if (opt == 15){
      //peanut
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
               double r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

               double theta = atan2(x[1],x[0]);
               
               double phi = acos(x[2]/r);

               if(r < 1e-9){
                S[i][j][k] = 0 - 0.5 - 0.2;
               }else{
                S[i][j][k] = r - 0.5 - 0.2 * sin(2*theta)*sin(phi); 
               }               
               
            }
         }
      }
   }
   else if (opt == 16){
      //popcorn
      double r0 = 0.6;
      double A = 2;
      double sigma = 0.2;
      array<array<double,3>,12> xk;

      for(int k = 0; k <=11; k++){
         if(k<=9){
            xk[k][0] = r0/sqrt(5.0) * 2 * cos(2 * k * M_PI/5 - floor(k/5.0) * M_PI);
            xk[k][1] = r0/sqrt(5.0) * 2 * sin(2 * k * M_PI/5 - floor(k/5.0) * M_PI);
            xk[k][2] = r0/sqrt(5.0) * pow(-1,floor(k/5.0));
         }else{
            xk[k][0] = r0 * 0;
            xk[k][1] = r0 * 0;
            xk[k][2] = r0 * pow(-1,k-10);
         }
      }
      for(int i = 0; i <= grid.nx[0]; i++){
         for(int j = 0; j <= grid.nx[1]; j++){
            for(int k = 0; k <= grid.nx[2]; k++){
               array<double,3> x = sub2coord(array<int,3>{i,j,k},grid);
               S[i][j][k] = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) - r0;
               double temp = 0.0;
               for(int n = 0; n <=11; n++){
                  temp += 2 * exp( -25 * pow(dist(x,xk[n]),2));
               }
               S[i][j][k] -= temp;
            }
         }
      }
   }else if (opt == 20){
      init_surf_protein_paper(S, grid);
   }
   else{
      cerr<<"invalid surface option";
      exit(1);
   }
}


void init_surf_perturb(double ***S, double radius, GridData &grid, int opt, PBData &pb, double tol){
  init_surf(S,  radius,  grid,  opt);
  if (globperturb == 0)
    perturb(S,grid.tol,pb,grid);
  else if (globperturb > 0)
    perturbstatus(S,grid.tol,globperturb,pb,grid);  
}

void cmdLine(int argc, char *argv[]){

//add option here
static struct option long_options[] = {
      {"GRIDNUM", required_argument, 0, 'n'},
      {"globcim", required_argument, 0, 'c'},
      {"globtestnum", required_argument, 0, 't'},
      {"SURFOPT", required_argument, 0, 's'},
      {"globtesta", required_argument, 0, 'a'},
      {"globlinsolve", required_argument, 0, 'l'},
      {"globcheck", required_argument, 0, 'k'},
      {"runStep", required_argument, 0, 'p'},
      {"globtime", required_argument, 0, 'i'},
      {"globheatsmooth", required_argument, 0, 'h'},
      {"eindex", required_argument, 0, 'x'},
      {"EPSILONM", required_argument, 0, 'M'},
      {"EPSILONP", required_argument, 0, 'P'},
      {"globperturb", required_argument, 0, 'r'},
      {"globwritemx", required_argument, 0, 'm'},
      {"globdist", required_argument, 0, 'd'},
      {"globsmall", required_argument, 0, 'S'},
      {"RADIUS", required_argument, 0, 'R'},
      {"globgridperturb", required_argument, 0, 'T'},
      {"tollinsolve", required_argument, 0, 'o'},
      {"PROTEINFILE", required_argument, 0, 'f'}
    };

   int c;
   int ei = 0;
   // add option here
   while ((c=getopt_long(argc,argv,"n:c:t:s:a:l:k:p:i:h:x:M:P:r:m:d:S:R:T:o:f:",long_options,NULL)) != -1){
      switch (c) {

         case 'n':// number of grid point
             GRIDNUM = atoi(optarg);
             break;

         case 'c':// cim method
             globcim = atoi(optarg);
             break;

         case 't':// test number 10 = litien's sin cos, 11 = cim paper
             globtestnum = atoi(optarg);
             break;

         case 's':// surface 0 circle, 1 litien's torus, 11 ellipsoid
             SURFOPT = atoi(optarg);
             break;

         case 'a':// 0 = without a term
             globtesta = atoi(optarg);
             break;

         case 'l':// linear solver
             globlinsolve = atoi(optarg);
             break;

         case 'k':// globcheck
             globcheck = atoi(optarg);
             break;

         case 'p':
             runStep = atoi(optarg);
             break;

         case 'i':
             globtime = atof(optarg);
             break;

         case 'h':
            globheatsmooth = atoi(optarg);
            break;

         case 'x':
            eindex[ei] = atoi(optarg);
            ei++;
            break;

         case 'M':
            EPSILONM = atoi(optarg);
         break;

         case 'P':
            EPSILONP = atoi(optarg);
         break;

         case 'r':
            globperturb = atoi(optarg);
         break;

         case 'm':
            globwritemx = atoi(optarg);
         break;

         case 'd':
            globdist = atoi(optarg);
         break;

         case 'S':
            globsmall = atoi(optarg);
         break;

         case 'R':
            RADIUS = atof(optarg);

          case 'T':
            globgridperturb = atoi(optarg);
         break;

         case 'o':
            tollinsolve = atof(optarg);
         break;

         case 'f':
            PROTEINFILE = string(optarg);
         break;
         
         default:
         // describe option here
             printf("Usage:\n"
               "-n --GRIDNUM\n"
               "-c --globcim\n"
               "-t --globtestnum: 10 = trig, 11 = cim paper>\n"
               "-s --SURFOPT: 0 = circle, 1 = wierd torus, 11 = ellipsoid, 12 = donut, 13 = banana, 14 = 8 balls, 15 = peanut, 16 = popcorn\n"
               "-a --globtesta, 0 = without\n"
               "-l --globlinsolve, 0 = BICGSTAB, 1 = ILU BICGSTAB, 2 = AMG, 3 = GMRES\n"
               "-k --globcheck, 1 = check convergence\n"
               "-p --runStep, number of step\n"
               "-i --globtime, tfinal\n"
               "-h --globheatsmooth, number heatsmooth step\n"
               "-x --eindex, index for debug\n"
               "-M --EPSILONM \n"
               "-P --EPSILONP \n"
               "-r --globperturb, perturb\n"
               "-m --globwritemx, write coupling matrix\n"
               "-d --globdist, differencing method\n"
               "-S --globsmall, store in small structure\n"
               "-R --radius <RADIUS>\n"
               "-T --globgridperturb\n"
               "-o --tollinsolve\n"
               "-f --PROTEINFILE\n");
               
             exit(-1);
      }
   }
   for (; optind < argc; optind++)
      printf ("Non-option argument %s\n", argv[optind]);
   
   if (globtestnum == 4){
      globtesta = 1;
   }

   if(globtesta == 1 && globcim == 7){
    cerr<<"globtesta=1 with globcim=7 not implemented"<<endl;
    exit(1);
   }
   // print option here
   printf("Problem: [globtestnum = %d]  [globtesta = %d] [EPSILONP = %f] [EPSILONM = %f] [PROTEINFILE = %s]\n", globtestnum, globtesta, EPSILONP, EPSILONM, PROTEINFILE.c_str());

   printf("Method: [globcim = %d] [globdist = %d] [globsmall = %d]  [globallcim345 = %d]\n", globcim,  globdist, globsmall,globallcim345);

   printf("Solver: [globlinsolve = %d] [tollinsolve = %e] \n",  globlinsolve, tollinsolve);

   printf("Geometry: [gridnum = %d] [SURFOPT = %d] [RADIUS = %f] [globperturb = %d] [globgridperturb = %d] \n", GRIDNUM, SURFOPT, RADIUS, globperturb, globgridperturb) ;

   printf("Motion: [globcheck = %d] [runStep = %d] [tfinal = %f] [globheatsmooth = %d]\n", globcheck, runStep, globtime, globheatsmooth);

   printf("Debug: [eindex %d %d %d] [globwritemx = %d] \n", eindex[0] , eindex[1], eindex[2], globwritemx);


   printf("Other options [globGSsmooth = %d] [globperturb = %d] [globdtdx2 = %d] [globheatsmoothtime = %f] [globbiasstat = %d]\n",
                        globGSsmooth,globperturb,globdtdx2,globheatsmoothtime,globbiasstat);
   
}

void checkwithexactvn(double ***vn, double ***S, PBData &pb, GridData &grid){
   double exactvn;
   if (globtestnum == 0){
      exactvn = 4.0*grid.radius/((1.0+grid.radius*grid.radius)*
                                 (1.0+grid.radius*grid.radius));
   }
   else if (globtestnum == 1){
      exactvn = 2.0*(1.0-pb.epsilonp/pb.epsilonm)*0.5*
                exp(2.0*(1.0-pb.epsilonp/pb.epsilonm)*grid.t);
   }
   else if (globtestnum == 2){
      exactvn = 2.0*grid.radius*exp(grid.radius*grid.radius)*
                (1.0-pb.epsilonp/pb.epsilonm);
   }
   else if (globtestnum == 3 || globtestnum == 4){
      exactvn = (1.0-pb.epsilonp/pb.epsilonm)*grid.radius/
                sqrt(grid.radius*grid.radius+fabs(1.0-pb.epsilonp/pb.epsilonm));
   }
   else{
      cerr<<"invalid motion test"<<endl;
      exit(1);
   }      


   vector<double> vnAll;
   vector<double> vnInterface;
   array<double,3> xAll; // coordinate of max vn extended
   array<double,3> xInterface; // coordinate of max vn interface
   double maxvnAll = 0; // max vn extended
   double maxvnInterface = 0; // max vn interface
   for(int i = 0; i <= grid.nx[0]; i++){
      for(int j = 0; j <= grid.nx[1]; j++){
         for(int k = 0; k <= grid.nx[2]; k++){
            array<int,3> x = {i, j, k};
            double absvn = abs(vn[i][j][k]);
            vnAll.push_back(absvn);
            if(absvn > maxvnAll){
               maxvnAll = absvn;
               xAll = sub2coord(array<int,3>{i,j,k},grid);
            }
            
            if(nearinterface(x, S, grid)){
               vnInterface.push_back(absvn);
               if(absvn > maxvnInterface){
                  maxvnInterface = absvn;
                  xInterface = sub2coord(array<int,3>{i,j,k},grid);
               }

            }
         }
      }
   }
   printf("abs(exactvn) = %f\n", abs(exactvn));
   printf("vn interface mean %f, variance %f, rmse %f\n", mean(vnAll), variance(vnAll), rmse(vnAll, abs(exactvn)));
   printf("maxvn interface %f at (%f,%f,%f)\n",maxvnInterface, xInterface[0], xInterface[1], xInterface[2]);
   printf("vn extend mean %f, variance %f, rmse %f\n", mean(vnInterface), variance(vnInterface), rmse(vnInterface, abs(exactvn)));
   printf("maxvn extend %f at (%f,%f,%f)\n",maxvnAll, xAll[0], xAll[1], xAll[2]);
   
}


// output ux at grid point and Du at interface
// modifed from getcim345Du, add *ux as argument
void getcim345DuAll(double &uint, double *ux, double *Du, int *index, int rstar, int sstar,
                 double ***u, double ***S, PBData &pb, GridData &grid)
{
   int r, s, t, i, m, n, sk, mid = 2, N = 2*mid;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];
   double d0[2*grid.dim];
   double ***dcoef[2*grid.dim];
   double **LU, **G;
   int gamma[grid.dim][2];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double alpha, beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;

   double D1u[grid.dim], ****D1ucoef, **D1uxcoef, **D1uxxcoef, 
          ***D1jumpuxxcoef;
   D1ucoef = new double ***[grid.dim];
   for (r = 0; r < grid.dim; r++)
      D1ucoef[r] = matrix(N,N,N);
   D1uxcoef = matrix(grid.dim-1,grid.dim-1);
   D1uxxcoef = matrix(grid.dim-1,grid.dim-1);
   D1jumpuxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);

   double finalalpha, finalD1u[grid.dim], ****finalD1ucoef, **finalD1uxcoef, 
          **finalD1uxxcoef;
   finalD1ucoef = new double ***[grid.dim];
   for (r = 0; r < grid.dim; r++)
      finalD1ucoef[r] = matrix(N,N,N);
   finalD1uxcoef = matrix(grid.dim-1,grid.dim-1);
   finalD1uxxcoef = matrix(grid.dim-1,grid.dim-1);

   double **D2u, *****D2ucoef, ***D2uxcoef, ***D2uxxcoef, **D2jumpuxxcoef;
   D2u = matrix(grid.dim-1,grid.dim-1);
   D2ucoef = new double ****[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      D2ucoef[r] = new double ***[grid.dim];
      for (s = 0; s < grid.dim; s++)
         D2ucoef[r][s] = matrix(N,N,N);
   }
   D2uxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   D2uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   D2jumpuxxcoef = matrix(grid.dim-1,grid.dim-1);

   double jumpD1u = 0.0, ***jumpD1ucoef, jumpD1uxcoef[grid.dim], 
          jumpD1uxxcoef[grid.dim], **jumpD1jumpuxxcoef;
   jumpD1ucoef = matrix(N,N,N);
   jumpD1jumpuxxcoef = matrix(grid.dim-1,grid.dim-1);

   double jumpD2u = 0.0, ***jumpD2ucoef, jumpD2uxcoef[grid.dim], 
          jumpD2uxxcoef[grid.dim];
   jumpD2ucoef = matrix(N,N,N);

   char yesD2[grid.dim][grid.dim];
   int eindex[grid.dim];
   eindex[0] = 17;
   eindex[1] = 5;
   eindex[2] = 19;

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
   {
      for (s = -1; s <= 1; s += 2)
      {
         rindex[r] = min(max(index[r]+s,0),grid.nx[r]);
         if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)
            gamma[r][(s+1)/2] = 1;
         else
            gamma[r][(s+1)/2] = 0;//gamma[r][s=0] = dim r, s=-1, has interface
      }
      rindex[r] = index[r];
   }

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   G = matrix(2*grid.dim-1,2*grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      dcoef[r] = matrix(N,N,N);

   if (evalarray(S,index) < 0.0)
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
   // get approximation of D2u at index
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         yesD2[m][n] = 0;
         if (!globdist)
            getcim345D2u(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],
                         D2jumpuxxcoef[m][n],yesD2[m][n],m,n,index,0,0,mid,S,grid);
         else
            getcim345D2udist(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],
                             D2jumpuxxcoef[m][n],yesD2[m][n],m,n,index,0,0,mid,S,grid);
      }

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
//            cout << "r = " << r << " " << sk << endl;
//            cout << "getting D2u" << endl;
            for (m = 0; m < grid.dim; m++)
               for (n = m+1; n < grid.dim; n++)
                  if (!globdist)
                     getcim345D2u(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                  D2uxxcoef[m][n],D2jumpuxxcoef[m][n],yesD2[m][n],
                                  m,n,index,r,sk,mid,S,grid);
                  else
                     getcim345D2udist(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                      D2uxxcoef[m][n],D2jumpuxxcoef[m][n],yesD2[m][n],
                                      m,n,index,r,sk,mid,S,grid);
            
            rindex[r] = index[r]+sk;
// getting interface info
            getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
            beta = 1.0-alpha;
// get derivatives
//            cout << "getting Du" << endl;
            getcim345Du(D1u,D1ucoef,D1uxcoef,D1uxxcoef,D1jumpuxxcoef,index,r,sk,alpha,
                        thesign,D2u,D2ucoef,D2uxcoef,D2uxxcoef,D2jumpuxxcoef,mid,grid);
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha,grid);
//            cout << "getting jump Du" << endl;
            getcim345jumpux(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                            jumpD1jumpuxxcoef,index,r,sk,alpha,thesign,normal,tangent,
                            mid,D1ucoef,D1uxcoef,D1uxxcoef,D1jumpuxxcoef,S,pb,grid);
//            cout << "getting jump D2u" << endl;
            getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha,thesign,normal,mid,D1u,D1ucoef,D1uxcoef,D1uxxcoef,
                             D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,D2uxxcoef,D2jumpuxxcoef,
                             jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                             jumpD1jumpuxxcoef,S,pb,grid);
            if (r == rstar && sk == sstar)
            {
               finalalpha = alpha;

               for (s = 0; s < grid.dim; s++)
                  finalD1u[s] = D1u[s];

               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 0;
               while (sindex[0] <= N)
               {
                  for (s = 0; s < grid.dim; s++)
                     setvalarray(finalD1ucoef[s],sindex,evalarray(D1ucoef[s],sindex));

                  (sindex[grid.dim-1])++;
                  for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
                  {
                     sindex[i] = 0;
                     (sindex[i-1])++;
                  }
               }
               for (s = 0; s < grid.dim; s++)
                  sindex[s] = mid;

               for (s = 0; s < grid.dim; s++)
                  for (t = 0; t < grid.dim; t++)
                  {
                     finalD1uxcoef[s][t] = D1uxcoef[s][t];
                     finalD1uxxcoef[s][t] = D1uxxcoef[s][t];
                  }
            }
// form d0 and dcoef's rth entry 
            d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                               sk*beta*jumpD1u/grid.dx[r]+
                                               0.5*beta*beta*jumpD2u)-
                                      sk*D1u[r]/grid.dx[r];

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                           thesign*(sk*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                    0.5*beta*beta*evalarray(jumpD2ucoef,sindex))-
                           sk*evalarray(D1ucoef[r],sindex)/grid.dx[r]);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)-
                        1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)+
                        1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid;
// form matrix's rth row
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m] = -thesign*
                                           (sk*beta*jumpD1uxxcoef[m]/grid.dx[r]+
                                            0.5*beta*beta*jumpD2uxxcoef[m])+
                                           sk*D1uxxcoef[r][m]/grid.dx[r];
            G[grid.dim*(sk+1)/2+r][r] += 0.5*(beta*beta-alpha*alpha);
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m+grid.dim] = -thesign*
                                                    (sk*beta*jumpD1uxcoef[m]/grid.dx[r]+
                                                     0.5*beta*beta*jumpD2uxcoef[m])+
                                                    sk*D1uxcoef[r][m]/grid.dx[r];
            rindex[r] = index[r];
         }
         else
         {
// if interior point
            d0[grid.dim*(sk+1)/2+r] = 0.0;
            for (s = 0; s < 2*grid.dim; s++)
               G[grid.dim*(sk+1)/2+r][s] = 0.0;
            G[grid.dim*(sk+1)/2+r][r] = 0.5;
            G[grid.dim*(sk+1)/2+r][grid.dim+r] = sk/grid.dx[r];
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,0.0);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,-1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid;
         }

// solve for uxx and put in matrix
   int j;
   double uxx[grid.dim];

   gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,2*grid.dim-1);
   for (j = 0; j < grid.dim; j++)
      uxx[j] = temp[j];
   for (j = 0; j < grid.dim; j++)
      ux[j] = temp[j+grid.dim];

   for (s = 0; s < grid.dim; s++)
      tindex[s] = index[s];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (n = 0; n < 2*grid.dim; n++)
         temp[n] = evalarray(dcoef[n],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-mid+sindex[s];
      for (j = 0; j < grid.dim; j++)
      {
         uxx[j] += temp[j]*evalarray(u,tindex);
         ux[j] += temp[j+grid.dim]*evalarray(u,tindex);
      }

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   if (globintorder == 3)
      uint = evalarray(u,index)+sstar*finalalpha*grid.dx[rstar]*ux[rstar]+
             0.5*(finalalpha*grid.dx[rstar])*(finalalpha*grid.dx[rstar])*uxx[rstar];
   else
      uint = evalarray(u,index)+sstar*finalalpha*grid.dx[rstar]*ux[rstar];

   for (s = 0; s < grid.dim; s++)
      Du[s] = finalD1u[s];

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-mid+sindex[s];
      for (s = 0; s < grid.dim; s++)
         Du[s] += evalarray(finalD1ucoef[s],sindex)*evalarray(u,tindex);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   for (s = 0; s < grid.dim; s++)
      for (t = 0; t < grid.dim; t++)
         Du[s] += finalD1uxcoef[s][t]*ux[t];

   for (s = 0; s < grid.dim; s++)
      for (t = 0; t < grid.dim; t++)
         Du[s] += finalD1uxxcoef[s][t]*uxx[t];


   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(G,2*grid.dim-1,2*grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(dcoef[r],N,N,N);

   for (r = 0; r < grid.dim; r++)
      free_matrix(D1ucoef[r],N,N,N);
   delete [] D1ucoef;
   free_matrix(D1uxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1jumpuxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);

   for (r = 0; r < grid.dim; r++)
      free_matrix(finalD1ucoef[r],N,N,N);
   delete [] finalD1ucoef;
   free_matrix(finalD1uxcoef,grid.dim-1,grid.dim-1);
   free_matrix(finalD1uxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(D2u,grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      for (s = 0; s < grid.dim; s++)
         free_matrix(D2ucoef[r][s],N,N,N);
      delete [] D2ucoef[r];
   }
   delete [] D2ucoef;
   free_matrix(D2uxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D2uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D2jumpuxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(jumpD1ucoef,N,N,N);
   free_matrix(jumpD1jumpuxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(jumpD2ucoef,N,N,N);
}

// check Du at gridpoint
void checkcim345DuAll(double ***u, double ***S, PBData &pb, GridData &grid)
{   
   double uint, ux[grid.dim], Du[grid.dim]; //results from getcim345
   
   array<int,3> maxErrIdx;
   array<int,3> maxErrItfIdx;
   double maxErr; // maxErr at all grid point
   double maxErrItf; // maxErr at grid point near interface
   vector<double> sqrdCompErrAll; // vector of error component at all grid point
   vector<double> sqrdCompErrItf; // vector of error component at grid point near interface
   
   double max_Dup_itf_err = 0.0;// error at interface plus side
   array<int,3> max_Dup_itf_err_idx;

   double max_Dup_itf_rel_err = 0.0;// relative error at interface plus side
   array<int,3> max_Dup_itf_rel_err_idx;
   
   double max_Dum_itf_err = 0.0;// error at interface minus side
   array<int,3> max_Dum_itf_err_idx;

   double max_Dum_itf_rel_err = 0.0;// relative error at interface minus side
   array<int,3> max_Dum_itf_rel_err_idx;
   
   // do not go through boundary
   for(int i = 1; i < grid.nx[0]; i++){
      for(int j = 1; j < grid.nx[1]; j++){
         for(int k = 1; k < grid.nx[2]; k++){
            int tindex[3] = {i,j,k};
            int rindex[3] = {i,j,k};

            int thestatus = getstatus5(S,tindex,grid);
            double thesign = (evalarray(S,tindex) < 0.0)? -1.0 : 1.0;
            
            double alpha, tangent[grid.dim], normal[grid.dim];
            if (thestatus>=2){ 
               //cim point, find rstar and sstar
               for (int r = 0; r < grid.dim; r++){
                  for (int s = -1; s <= 1; s += 2){
                     rindex[r] = tindex[r]+s;

                     if ((evalarray(S,tindex) < 0.0)+(evalarray(S,rindex) < 0.0) == 1){
                        
                        getcim345DuAll(uint,ux,Du,tindex,r,s,u,S,pb,grid);
                        
                        getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
                        
                        //Du exact at interface
                        array<double,3> du_itf_exact {getDu(tindex,0,r,s,alpha,thesign,grid),
                                                getDu(tindex,1,r,s,alpha,thesign,grid),
                                                getDu(tindex,2,r,s,alpha,thesign,grid)};
                        double max_Du_itf_exact = max(max(fabs(du_itf_exact[0]),fabs(du_itf_exact[1])),fabs(du_itf_exact[2]));
                        //Du at interface

                        array<double,3> Du_itf_err =  { fabs(Du[0] - du_itf_exact[0]),
                                                   fabs(Du[1] - du_itf_exact[1]),
                                                   fabs(Du[2] - du_itf_exact[2])};
                        double max_Du_itf_err = max(max(Du_itf_err[0],Du_itf_err[1]),Du_itf_err[2]);
                        
                        double Du_itf_rel_err = max_Du_itf_err/max_Du_itf_exact; //relative error
                        
                        if(thesign<0){
                           max_Dum_itf_err = (max_Du_itf_err > max_Dum_itf_err)? max_Du_itf_err : max_Dum_itf_err;

                           max_Dum_itf_err_idx = array<int,3>{i,j,k};

                           max_Dum_itf_rel_err = (Du_itf_rel_err > max_Dum_itf_rel_err)? Du_itf_rel_err : max_Dum_itf_rel_err;
                           max_Dum_itf_rel_err_idx = array<int,3>{i,j,k};
                        }else{
                           max_Dup_itf_err = (max_Du_itf_err > max_Dup_itf_err)? max_Du_itf_err : max_Dup_itf_err;
                           max_Dup_itf_err_idx = array<int,3>{i,j,k};

                           max_Dup_itf_rel_err = (Du_itf_rel_err > max_Dup_itf_rel_err)? Du_itf_rel_err : max_Dup_itf_rel_err;
                           max_Dup_itf_rel_err_idx = array<int,3>{i,j,k};
                        }
                     }else{
                        continue;
                     }
                  }//end of s loop
                  rindex[r] = tindex[r];
               }// end of r loop


            }else if (thestatus == 1){ 
               //interior, central differencing
               ux[0] = (u[i+1][j][k] - u[i-1][j][k]) / (2*grid.dx[0]);
               ux[1] = (u[i][j+1][k] - u[i][j-1][k]) / (2*grid.dx[1]);
               ux[2] = (u[i][j][k+1] - u[i][j][k-1]) / (2*grid.dx[2]);

            }else{
               cerr<<"check Du unclassified"<<endl;
               exit(1);
            }

            

            double xerr = abs(ux[0] - getDu(tindex,0, 0, 0, 0.0, thesign, grid));
            double yerr = abs(ux[1] - getDu(tindex,1, 0, 0, 0.0, thesign, grid));
            double zerr = abs(ux[2] - getDu(tindex,2, 0, 0, 0.0, thesign, grid));
            double tmpMaxErr = max(max(xerr,yerr),zerr);

            if (tmpMaxErr>maxErr){
               maxErr = tmpMaxErr;
               maxErrIdx = array<int,3>{i,j,k};
            }
            sqrdCompErrAll.push_back(xerr);
            sqrdCompErrAll.push_back(yerr);
            sqrdCompErrAll.push_back(zerr);

            if(thestatus>=2){
               if (tmpMaxErr>maxErrItf){
                  maxErrItf = tmpMaxErr;
                  maxErrItfIdx = array<int,3>{i,j,k};
               }
               sqrdCompErrItf.push_back(xerr);
               sqrdCompErrItf.push_back(yerr);
               sqrdCompErrItf.push_back(zerr);
               
            }
         }//end k
      }// end j
   }// end i
    
   cout<<"[checkcim345DuAll]"<<endl;
   printf("max Du err All = %f at (%d,%d,%d)\n",maxErr,maxErrIdx[0],maxErrIdx[1],maxErrIdx[2]);
   printf("max Du err Interface grid point = %f at (%d,%d,%d)\n",maxErrItf,maxErrItfIdx[0],maxErrItfIdx[1],maxErrItfIdx[2]);
   printf("rmse Du err All = %f \n",rmse(sqrdCompErrAll));
   printf("rmse Du err Interface grid point = %f \n",rmse(sqrdCompErrItf));

   printf("max Du err negative interface = %f at (%d,%d,%d)\n", max_Dum_itf_err, max_Dum_itf_err_idx[0], max_Dum_itf_err_idx[1], max_Dum_itf_err_idx[2]);
   printf("max Du err positive interface = %f at (%d,%d,%d)\n", max_Dup_itf_err, max_Dup_itf_err_idx[0], max_Dup_itf_err_idx[1], max_Dup_itf_err_idx[2]);

   printf("max Du rel err negative interface = %f at (%d,%d,%d)\n", max_Dum_itf_rel_err, max_Dum_itf_rel_err_idx[0], max_Dum_itf_rel_err_idx[1], max_Dum_itf_rel_err_idx[2]);
   printf("max Du rel err positive interface = %f at (%d,%d,%d)\n", max_Dup_itf_rel_err, max_Dup_itf_rel_err_idx[0], max_Dup_itf_rel_err_idx[1], max_Dup_itf_rel_err_idx[2]);
}



void read_protein(double ***S, GridData &grid, vector< vector<double> >& p){
  cout<<"using protein file "<<PROTEINFILE<<endl;
   ifstream infile( PROJECTDIR + PROTEINFILE);
   if(!infile.is_open()){
      cerr<<"file do not exist"<<endl;
      exit(1);
   }

   string line;
   while(getline(infile,line)){
      stringstream ss(line);
      string x,y,z,r;
      getline(ss,x,',');
      getline(ss,y,',');
      getline(ss,z,',');
      getline(ss,r,',');
   
      
      p.push_back(std::vector<double> {stod(x),stod(y),stod(z),stod(r)});
   }

   for(int i = 0; i < 3; i++){
      printf("%d %f %f %f %f\n",i,p[i][0],p[i][1],p[i][2],p[i][3]);
   }

   for(int i = p.size()-3; i < p.size(); i++){
      printf("%d %f %f %f %f\n",i,p[i][0],p[i][1],p[i][2],p[i][3]);
   }   
}

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

inline double chi(double x){
   const double eta = 1.0/ETA; //in paper 1/40
   return 1.0/2.0 * (1.0 + tanh(x/eta));
}

void init_surf_protein_paper(double ***S, GridData& grid){
  // if file exist
  string file = PROTEINDIR + string("surf_protein_paper") + "_eta"+to_string( (int)ETA) + "_gd" + to_string(GRIDNUM) + "_" + PROTEINFILE;
  if(fileExist(file)){
    cout<<"reading "<< file<<endl;
    readMatrix(file, S, grid);
    return;
  }
  cout<<"file "<<file <<" do not exist, construct from protein"<<endl;
  // if file do not exist
  vector< vector<double> > p;
  read_protein(S,grid,p);
  const double c = 0.25;
   #pragma omp parallel for collapse(3)
   for(int i = 0; i <= grid.nx[0]; i++){
      for(int j = 0; j <= grid.nx[1]; j++){
         for(int k = 0; k <= grid.nx[2]; k++){
            array<double,3> x = sub2coord(array<int,3> {i,j,k}, grid);
            double sum = 0.0;

            for(int n = 0; n < p.size(); n++){
               array<double,3> px = {p[n][0],p[n][1],p[n][2]};//px is coordinate
               sum += chi(p[n][3] - dist(x,px)); //p[n][3] = radius of atom n
            }
            S[i][j][k] = c - sum;
         }
      }
   }

   // FillVoid(S, grid);
   // write to file
   cout<<"write "<< file<<endl;
   write_field(file, S, grid);

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

// get status of point index: 0 = boundary, 1 = interior, 2 = cim3, 3 = cim5, 4 = cim4, 5 = cim1
int getstatus5debug(double ***S, int *index, GridData &grid)
{
   int r, sk, m, n, mid = 1, rindex[grid.dim], sk2[4];
   int yesinterior, thecim, cimstatus[grid.dim][grid.dim];
   double ***junk1, junk2[grid.dim], junk3[grid.dim];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
      {
         rindex[r] = index[r]+sk;
         if (rindex[r] < 0 || rindex[r] > grid.nx[r])
            return 0;//if index is boundary
         rindex[r] = index[r];
      }

   yesinterior = 1;
   for (r = 0; r < grid.dim; r++)
   {
      for (sk = -1; sk <= 1; sk += 2)
      {
         rindex[r] = index[r]+sk;
         if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1) // if near interface
            yesinterior = 0;
      }
      rindex[r] = index[r];
   }
   if (yesinterior)
      return 1;

// old incorrect version returned cim 5 if cim 5 for any m and n 
   junk1 = matrix(2,2,2);
//   thecim = 3;
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         cimstatus[m][n] = 3; 
         if (!yessk2(sk2,m,n,index,S,grid))//if do not have usual mixed derivative
         {

            if (yescim5D2(junk1,junk2,junk3,m,n,index,mid,S,grid))//if has cim5 mixed derivative
//               if (thecim == 3)
//                  thecim = 5;
               cimstatus[m][n] = 5; 
            else 
            {
//               if (thecim == 3 || thecim == 5)
//                  thecim = 4;
               cimstatus[m][n] = 4; 
               for (r = 0; r < grid.dim; r++)
                  for (sk = -1; sk <= 1; sk += 2)
                  {
                     rindex[r] = index[r]+sk;
                     if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1 &&
                         !yessk2(sk2,m,n,rindex,S,grid))
//                        thecim = 0;
                        cout<<"m = "<<m<<", n = "<<n<<", r = "<<r<<", sk = "<<sk<<", cimstatus = "<<cimstatus[m][n]<<endl;
                        cimstatus[m][n] = 0;
                     rindex[r] = index[r];
                  }
            }
         }
      }
   free_matrix(junk1,2,2,2);

   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++){
         cout<<"m = "<<m<<", n = "<<n<<", cimstatus = "<<cimstatus[m][n]<<endl;
      }



   thecim = 3;
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
         if (cimstatus[m][n] == 0)
            return 5;
         else if (cimstatus[m][n] == 4)
            thecim = 4;
         else if (thecim == 3 && cimstatus[m][n] == 5)
            thecim = 5;
   
   cout<<"thecim = "<<thecim<<endl;

   if (thecim == 3)
   {
// second derivatives exhibit cim 3 structure but have to still check first derivatives
      int count = 0;

      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];
      for (r = 0; r < grid.dim && count < 2; r++)
      {
         count = 0;
         for (sk = -1; sk <= 1; sk += 2)
         {
            rindex[r] = index[r]+sk;
            if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)//if change sign, count++
               count++;
         }
         rindex[r] = index[r];
      }
      if (count < 2)
         return 2;
      else
         thecim = 5;
   }

   if (thecim == 5)
      return 3;
   else if (thecim == 4)
      return 4;
   else
      return 5; 
}


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

// estimate condition number
double condest(const vector<vector<double>> &A){
   int n = A.size();
   double n1 = norm1(A);
   vector<double> x(n, 1.0/n);
   int i1 = -1, i2;
   double c1 = 0, c2 = 0;
   

   while(true){
      x = gesolve(A,x);

      c2 = sum(abs(x));

      x = sign(x);

      x = gesolve(transpose(A),x);

      i2 = max_abs_idx(x);

      if(0 <= i1){
         if(i1==i2 || c2<=c1){
            break;
         }
      }

      i1 = i2;
      c1 = c2;

      fill(x.begin(),x.end(),0);
      x[i1] = 1.0;
   }

   return c2*n1;
}

double condest(double **G, int n){
   return condest(Pointer2dToVector(G,n,n));
}

// to use without a(x), pass ***a as nullptr
// only difference is getcim345jumpuxx
void gmatrix(double **G, int *index, int gamma[][2], double ***S, double***a, PBData &pb, GridData &grid){
   //copy every variable excetp G
   int r, s, t, i, j, m, n, sk, mid = 2, N = 2*mid;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim], Narray[grid.dim];
   double d0[2*grid.dim];
   double ***dcoef[2*grid.dim];
   double **LU, value;
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double **alpha = matrix(grid.dim-1,1), beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;

            double exactd[2*grid.dim], exactres, tempsign;
            int tempindex[grid.dim];

   for (r = 0; r < grid.dim; r++)
      Narray[r] = N;

   double ***D1u, ******D1ucoef, ****D1uxcoef, ****D1uxxcoef, ***D1jumpuxxcoef;
   D1u = new double **[grid.dim];
   D1ucoef = new double *****[grid.dim];
   D1uxcoef = new double ***[grid.dim];
   D1uxxcoef = new double ***[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      D1u[r] = new double *[2];
      D1ucoef[r] = new double ****[2];
      D1uxcoef[r] = new double **[2];
      D1uxxcoef[r] = new double **[2];
      for (s = 0; s <= 1; s++)
      {
         D1u[r][s] = NULL;
         D1ucoef[r][s] = NULL;
         D1uxcoef[r][s] = NULL;
         D1uxxcoef[r][s] = NULL;
      }
   }
   D1jumpuxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   
   double **D2u, *****D2ucoef, ***D2uxcoef, ***D2uxxcoef, **D2jumpuxxcoef;
   D2u = matrix(grid.dim-1,grid.dim-1);
   D2ucoef = new double ****[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      D2ucoef[r] = new double ***[grid.dim];
      for (s = 0; s < grid.dim; s++)
         D2ucoef[r][s] = matrix(N,N,N);
   }
   D2uxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   D2uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   D2jumpuxxcoef = matrix(grid.dim-1,grid.dim-1);

   double jumpD1u = 0.0, ***jumpD1ucoef, jumpD1uxcoef[grid.dim], 
          jumpD1uxxcoef[grid.dim], **jumpD1jumpuxxcoef;
   jumpD1ucoef = matrix(N,N,N);
   jumpD1jumpuxxcoef = matrix(grid.dim-1,grid.dim-1);

   double jumpD2u = 0.0, ***jumpD2ucoef, jumpD2uxcoef[grid.dim], 
          jumpD2uxxcoef[grid.dim];
   jumpD2ucoef = matrix(N,N,N);

   double ux[grid.dim], ****uxcoef;
   uxcoef = new double ***[grid.dim];
   for (r = 0; r < grid.dim; r++)
      uxcoef[r] = matrix(N,N,N);

   char yesD2[grid.dim][grid.dim];
   
   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
   {
      cout << "S" << endl;
      for (t = -1; t <= 1; t++) 
      {
         for (s = 1; s >= -1; s--) 
         {
            for (r = -1; r <= 1; r++) 
               printf("%4.16f ",S[index[0]+r][index[1]+s][index[2]+t]);
            cout << endl;
         }
         cout << endl;
      }
      cout << "u" << endl;
      for (t = -1; t <= 1; t++) 
      {
         for (s = 1; s >= -1; s--) 
         {
            for (r = -1; r <= 1; r++) 
            {
               tindex[0] = index[0]+r;
               tindex[1] = index[1]+s;
               tindex[2] = index[2]+t;
               thesign = (evalarray(S,index) < 0.0)?-1:1;
               printf("%4.16f ",getu(tindex,0,0,0.0,thesign,grid));
            }
            cout << endl;
         }
         cout << endl;
      }
   }

   value = 0.0;

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      dcoef[r] = matrix(N,N,N);


   // actual calculation
   if (evalarray(S,index) < 0.0)
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

   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
// yesD2 only defined for m < n
         yesD2[m][n] = 0;
         if (!globdist)
            getcim345D2u(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],
                         D2jumpuxxcoef[m][n],yesD2[m][n],m,n,index,0,0,mid,S,grid);
         else
            getcim345D2udist(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],
                             D2jumpuxxcoef[m][n],yesD2[m][n],m,n,index,0,0,mid,S,grid);
      }

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
            D1u[r][(sk+1)/2] = new double[grid.dim];
            D1ucoef[r][(sk+1)/2] = new double ***[grid.dim];
            for (t = 0; t < grid.dim; t++)
               D1ucoef[r][(sk+1)/2][t] = matrix(N,N,N);
            D1uxcoef[r][(sk+1)/2] = matrix(grid.dim-1,grid.dim-1);
            D1uxxcoef[r][(sk+1)/2] = matrix(grid.dim-1,grid.dim-1);

            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               cout<<endl<<"dim = "<<r <<" sk = "<<sk << endl;
            }

            for (m = 0; m < grid.dim; m++)
               for (n = m+1; n < grid.dim; n++)
               {
                  if (!globdist)
                     getcim345D2u(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                  D2uxxcoef[m][n],D2jumpuxxcoef[m][n],yesD2[m][n],
                                  m,n,index,r,sk,mid,S,grid);
                  else
                     getcim345D2udist(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                      D2uxxcoef[m][n],D2jumpuxxcoef[m][n],yesD2[m][n],
                                      m,n,index,r,sk,mid,S,grid);
                  if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
                  {
                     cout <<"computed D2u in ("<< m<<","<<n<<") plane" << endl;
                     cout << " apprx = "
                          << evalcoef(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                      D2uxxcoef[m][n],index,0,0,0.0,mid,thesign,grid) 
                          << " exact = "
                          << getD2u(index,m,n,r,sk,alpha[r][(sk+1)/2],thesign,grid)
                          << " error = "
                          << evalcoef(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                      D2uxxcoef[m][n],index,0,0,0.0,mid,thesign,grid)-
                             getD2u(index,m,n,r,sk,alpha[r][(sk+1)/2],thesign,grid) 
                          << endl;
                  }
               }
            
            rindex[r] = index[r]+sk;
// getting interface info
            getinterfaceinfo(alpha[r][(sk+1)/2],tangent,normal,S,index,rindex,grid);
            beta = 1.0-alpha[r][(sk+1)/2];
// get derivatives
//            cout << "getting Du" << endl;
            getcim345Du(D1u[r][(sk+1)/2],D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                        D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,index,r,sk,
                        alpha[r][(sk+1)/2],thesign,D2u,D2ucoef,D2uxcoef,D2uxxcoef,
                        D2jumpuxxcoef,mid,grid);
            double rhs = 0.0;
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               cout << "computed Du" << endl;
               for (m = 0; m < grid.dim; m++)
                  cout<<"dim "<< m <<" apprx = " 
                      << evalcoef(D1u[r][(sk+1)/2][m],D1ucoef[r][(sk+1)/2][m],
                                   D1uxcoef[r][(sk+1)/2][m],D1uxxcoef[r][(sk+1)/2][m],
                                   index,0,0,0.0,mid,thesign,grid) 
                      << ", exact ="
                      << getDu(index,m,r,sk,alpha[r][(sk+1)/2],thesign,grid)
                      << ", error ="
                      << evalcoef(D1u[r][(sk+1)/2][m],D1ucoef[r][(sk+1)/2][m],
                                   D1uxcoef[r][(sk+1)/2][m],D1uxxcoef[r][(sk+1)/2][m],
                                   index,0,0,0.0,mid,thesign,grid)-
                          getDu(index,m,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
               
               rhs += sk*grid.dx[r]*
                      evalcoef(D1u[r][(sk+1)/2][r],D1ucoef[r][(sk+1)/2][r],
                               D1uxcoef[r][(sk+1)/2][r],D1uxxcoef[r][(sk+1)/2][r],
                               index,0,0,0.0,mid,thesign,grid);
               cout << "rhs 1 = " << rhs << ", exact = "
                    << sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) //rhs 1 = h ux-
                    << endl;
            }
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha[r][(sk+1)/2],grid);
//            cout << "getting jump Du" << endl;
            getcim345jumpux(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                            jumpD1jumpuxxcoef,index,r,sk,alpha[r][(sk+1)/2],thesign,
                            normal,tangent,mid,D1ucoef[r][(sk+1)/2],
                            D1uxcoef[r][(sk+1)/2],D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,
                            S,pb,grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               double x[grid.dim];
               sub2coord(x,index,grid);
               x[r] += sk*alpha[r][(sk+1)/2]*grid.dx[r];
               cout << "computed jump in D1u = "
                    << evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef, 
                                jumpD1jumpuxxcoef,index,0,0,0.0,S,grid) << " "
                    <<", exact = "
                    << getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                       getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid) 
                    // << "another exact "
                    // << -(pb.epsilonp-pb.epsilonm)/ethere*
                    //     (getDu(index,0,r,sk,alpha[r][(sk+1)/2],thesign,grid)*normal[0]+
                    //      getDu(index,1,r,sk,alpha[r][(sk+1)/2],thesign,grid)*normal[1]+
                    //      getDu(index,2,r,sk,alpha[r][(sk+1)/2],thesign,grid)*normal[2])*
                    //     normal[r] << " "
                    // << 2.0*x[r]*(1.0-pb.epsilonp/pb.epsilonm) 
                    << " error =  "
                    << evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                                jumpD1jumpuxxcoef,index,0,0,0.0,S,grid)-
                       (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid)) << endl;

               rhs += -thesign*sk*beta*grid.dx[r]*
                      evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                               jumpD1jumpuxxcoef,index,0,0,0.0,S,grid);
               cout << "rhs 2 = " << rhs << " error =  "
                    << sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid)) << endl; //rhs2 = h ux- + beta h [ux]
            }
//            cout << "getting jump D2u" << endl;
            
            if(a){
              // if a is not nullptr
              getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha[r][(sk+1)/2],thesign,normal,mid,a,D1u[r][(sk+1)/2],
                             D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                             D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
                             D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                             jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);
            }else{
              getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha[r][(sk+1)/2],thesign,normal,mid,D1u[r][(sk+1)/2],
                             D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                             D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
                             D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                             jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);  
            }
            
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               cout << "computed jump in D2u = "
                    << evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,0,0,
                                0.0,mid,thesign,grid) 
                    << ", exact = "
                    << getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid)
                    << ", error = "
                    << evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,0,0,
                                0.0,mid,thesign,grid)-
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid)) 
                    << endl;
               rhs += -thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                      evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,0,0,
                               0.0,mid,thesign,grid); 
               cout << "rhs 3 = " << rhs << " " //rhs4 = h ux- + beta h [ux] + (1/2) (betta h)^2  [uxx]
                    << ", exact = "
                    << -thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid)) << endl;
               
               rhs += 0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2])*grid.dx[r]*grid.dx[r]*
                      getD2u(index,r,r,0,0,0.0,thesign,grid); // uxx in dim r
               
               cout << "rhs 4 = " << rhs << " " //rhs4 = h ux- + beta h [ux] + (1/2) (beta h)^2  [uxx] + (1/2) h^2 (beta^2-alpha^2) (uxx-)
                    << ", exact = "
                    << -thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2])*grid.dx[r]*grid.dx[r]*
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
               
               cout << 0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2])*grid.dx[r]*grid.dx[r]*getD2u(index,r,r,0,0,0.0,thesign,grid) << " "
                    << 0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2])*grid.dx[r]*grid.dx[r]*getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << " "
                    << getD2u(index,r,r,0,0,0.0,thesign,grid) << " "
                    << getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
               
               int tempindex[grid.dim];
               for (m = 0; m < grid.dim; m++)
                  tempindex[m] = index[m];
               tempindex[r] += sk; //tempindex at the other side of interface
               cout << "rhs = " << rhs << endl;
                      // asume cim12 p25, (u-) - alpha (ux-) + 0.5 (alpha h)^2 uxx-
               cout << getu(index,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       sk*alpha[r][(sk+1)/2]*grid.dx[r]*
                       getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)+
                       0.5*alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2]*grid.dx[r]*grid.dx[r]*
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       
                       (getu(tempindex,r,-sk,beta,-thesign,grid)+
                        sk*beta*grid.dx[r]*
                        getDu(tempindex,r,r,-sk,beta,-thesign,grid)+
                        0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                        getD2u(tempindex,r,r,r,-sk,beta,-thesign,grid));
               cout << ", should ="
                    << thesign*tau
                       -thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       0.5*grid.dx[r]*grid.dx[r]*(beta*beta-alpha[r][(sk+1)/2]*
                                                            alpha[r][(sk+1)/2])*
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
              
               cout << "thesign = " << thesign << " and dim = " << r << " and direction = " << sk << endl;
               cout << " move to same side should be zero "//move everything to one side
                    << getu(index,0,0,0.0,thesign,grid)-
                       getu(index,r,sk,1.0,-thesign,grid)
                       + (-thesign)*tau
                       -thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       0.5*grid.dx[r]*grid.dx[r]*(beta*beta-alpha[r][(sk+1)/2]*
                                                            alpha[r][(sk+1)/2])*
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
               cout << "(at grid) u_here - u_there "
                    << " = " << getu(index,0,0,0.0,thesign,grid) 
                    << " - " << getu(index,r,sk,1.0,-thesign,grid)
                    << " = "<<getu(index,0,0,0.0,thesign,grid) - getu(index,r,sk,1.0,-thesign,grid) << endl;
               cout << "u jump = "
                    << getu(index,r,sk,alpha[r][(sk+1)/2],1,grid) - getu(index,r,sk,alpha[r][(sk+1)/2],-1,grid) << endl;
               double x[grid.dim];
               sub2coord(x,index,grid);
               x[r] += sk*alpha[r][(sk+1)/2]*grid.dx[r];
               cout << "error in radius = " <<sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-grid.radius0 << endl;
            }
// form d0 and dcoef's rth entry, d0 + dcoeff u = G [ux uxx]
            d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                               sk*beta*jumpD1u/grid.dx[r]+
                                               0.5*beta*beta*jumpD2u)-
                                      sk*D1u[r][(sk+1)/2][r]/grid.dx[r];
            exactd[grid.dim*(sk+1)/2+r] = d0[grid.dim*(sk+1)/2+r];//constant part is exact

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            { 
                //coefficient of u[sindex], move rhs to lhs
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                           thesign*(sk*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                    0.5*beta*beta*evalarray(jumpD2ucoef,sindex))-
                           sk*evalarray(D1ucoef[r][(sk+1)/2][r],sindex)/grid.dx[r]);
               for (m = 0; m < grid.dim; m++)
                  tempindex[m] = index[m]+sindex[m]-mid;
               if (evalarray(S,tempindex) < 0.0)
                  tempsign = -1.0;
               else
                  tempsign = 1.0;
               exactd[grid.dim*(sk+1)/2+r] += evalarray(dcoef[grid.dim*(sk+1)/2+r],
                                                        sindex)*
                                              getu(tempindex,0,0,0.0,tempsign,grid);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
             // coeff due to ( u_{i,j+1} - u_{i,j} )/ h^2 
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)-
                        1.0/(grid.dx[r]*grid.dx[r]));
            for (m = 0; m < grid.dim; m++)
               tempindex[m] = index[m]+sindex[m]-mid;
            if (evalarray(S,tempindex) < 0.0)
               tempsign = -1.0;
            else
               tempsign = 1.0;
            exactd[grid.dim*(sk+1)/2+r] += -1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(tempindex,0,0,0.0,tempsign,grid);
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)+
                        1.0/(grid.dx[r]*grid.dx[r]));
            for (m = 0; m < grid.dim; m++)
               tempindex[m] = index[m]+sindex[m]-mid;
            if (evalarray(S,tempindex) < 0.0)
               tempsign = -1.0;
            else
               tempsign = 1.0;
            exactd[grid.dim*(sk+1)/2+r] += 1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(tempindex,0,0,0.0,tempsign,grid);
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               double somezeros[grid.dim];
               for (m = 0; m < grid.dim; m++)
                  somezeros[m] = 0.0;
               cout << "Row = " <<grid.dim*(sk+1)/2+r
                    << ", exactd = "
                    << exactd[grid.dim*(sk+1)/2+r]
                    << ", approx = "
                    << evalcoef(d0[grid.dim*(sk+1)/2+r],dcoef[grid.dim*(sk+1)/2+r],
                                somezeros,somezeros,index,0,0,0.0,mid,S,grid)
                    << endl;
            }
            sindex[r] = mid;
// form matrix's rth row
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m] = -thesign*
                                           (sk*beta*jumpD1uxxcoef[m]/grid.dx[r]+
                                            0.5*beta*beta*jumpD2uxxcoef[m])+
                                           sk*D1uxxcoef[r][(sk+1)/2][r][m]/grid.dx[r]; // coeff of u_{xx} from [u_xx] [u_x]

            G[grid.dim*(sk+1)/2+r][r] += 0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2]); // 1/2 (beta^2-alpha^2)u_{xx}
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m+grid.dim] = -thesign*
                                                    (sk*beta*jumpD1uxcoef[m]/grid.dx[r]+
                                                     0.5*beta*beta*jumpD2uxcoef[m])+
                                                    sk*D1uxcoef[r][(sk+1)/2][r][m]/
                                                    grid.dx[r];// coeff of u_x from [u_xx] [u_x]
            if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
            {
               exactres = exactd[grid.dim*(sk+1)/2+r];
               for (m = 0; m < 2*grid.dim; m++)
                  if (m < grid.dim)
                     exactres -= G[grid.dim*(sk+1)/2+r][m]*
                                 getD2u(index,m,m,0,0,0.0,thesign,grid);
                  else
                     exactres -= G[grid.dim*(sk+1)/2+r][m]*
                                 getDu(index,m-grid.dim,0,0,0.0,thesign,grid);
               cout << "d_exact - G_apprx DDu_exact = " << exactres << endl;
               cout << "tau = "
                    << getu(index,r,sk,alpha[r][(sk+1)/2],1,grid)-
                       getu(index,r,sk,alpha[r][(sk+1)/2],-1,grid)
                    << ", tau/h^2  ="
                    << (getu(index,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getu(index,r,sk,alpha[r][(sk+1)/2],-1,grid))/
                       (grid.dx[r]*grid.dx[r]) << endl;
               cout << "-ehere lap(u) - f = "
                    << -ehere*(getD2u(index,0,0,0,0,0.0,thesign,grid)+
                               getD2u(index,1,1,0,0,0.0,thesign,grid)+
                               getD2u(index,2,2,0,0,0.0,thesign,grid))-
                       getf(index,0,0,0.0,thesign,pb,grid) 
                    << " =  "
                    << -ehere*(getD2u(index,0,0,0,0,0.0,thesign,grid)+
                               getD2u(index,1,1,0,0,0.0,thesign,grid)+
                               getD2u(index,2,2,0,0,0.0,thesign,grid)) 
                    << " - "
                    << getf(index,0,0,0.0,thesign,pb,grid) 
                    << endl
                    << "uxx = "<<getD2u(index,0,0,0,0,0.0,thesign,grid) << " "
                    << "uyy = "<<getD2u(index,1,1,0,0,0.0,thesign,grid) << " "
                    << "uzz = "<<getD2u(index,2,2,0,0,0.0,thesign,grid) << " "
                    << "ehere = "<<ehere 
                    << endl;
               double x[grid.dim];
               sub2coord(x,index,grid);
               x[r] += sk*alpha[r][(sk+1)/2]*grid.dx[r];
               cout <<"at interface x = ("<< x[0] << ", " << x[1] << ", " << x[2] <<")"<< endl;
               cout <<"error of radius" << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-grid.radius
                    << ", alpha = "<< alpha[r][(sk+1)/2] << endl;
            }
            rindex[r] = index[r];
         }
         else
         {
// if interior point
            d0[grid.dim*(sk+1)/2+r] = 0.0;
            exactd[grid.dim*(sk+1)/2+r] = d0[grid.dim*(sk+1)/2+r];
            for (s = 0; s < 2*grid.dim; s++)
               G[grid.dim*(sk+1)/2+r][s] = 0.0;
            G[grid.dim*(sk+1)/2+r][r] = 0.5;
            G[grid.dim*(sk+1)/2+r][grid.dim+r] = sk/grid.dx[r];
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,0.0);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,-1.0/(grid.dx[r]*grid.dx[r]));
            exactd[grid.dim*(sk+1)/2+r] += -1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(index,0,0,0.0,thesign,grid);
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
            exactd[grid.dim*(sk+1)/2+r] += 1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(index,r,sk,1.0,thesign,grid);
            sindex[r] = mid;
         }

   double uxval[grid.dim], uxxval[grid.dim]; // for chekcing Dusmall

   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2])
   {
      cout<<endl<<"G matrix"<<endl;
      for (m = 0; m < 2*grid.dim; m++)
      {
         for (n = 0; n < 2*grid.dim; n++)
            cout << setw(10)<< G[m][n] << " ";
         cout << endl;
      }
      cout<<"exactd"<<endl;
      for (m = 0; m < 2*grid.dim; m++)
         cout << exactd[m] << endl;
      gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
      forwardbacksub0(exactd,exactd,LU,PLR,PLC,2*grid.dim-1);

      cout<<"solved | exact"<<endl;
      for (m = 0; m < 2*grid.dim; m++)
      {
         cout  << exactd[m] << " | ";
         if (m < grid.dim)
            cout << getD2u(index,m,m,0,0,0.0,thesign,grid) << endl;
         else
            cout << getDu(index,m-grid.dim,0,0,0.0,thesign,grid) << endl;
      }
      cout <<" solved -eps lap(u) - f =" <<-ehere*(exactd[0]+exactd[1]+exactd[2])-getf(index,0,0,0.0,thesign,pb,grid)
           << " =  "<< -ehere*(exactd[0]+exactd[1]+exactd[2]) 
           << " - " << getf(index,0,0,0.0,thesign,pb,grid) << endl;
      
   }

   free_matrix(alpha,grid.dim-1,1);
   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(dcoef[r],N,N,N);

   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
            delete [] D1u[r][(sk+1)/2];
            for (t = 0; t < grid.dim; t++)
               free_matrix(D1ucoef[r][(sk+1)/2][t],N,N,N);
            delete [] D1ucoef[r][(sk+1)/2];
            free_matrix(D1uxcoef[r][(sk+1)/2],grid.dim-1,grid.dim-1);
            free_matrix(D1uxxcoef[r][(sk+1)/2],grid.dim-1,grid.dim-1);
         }
   for (r = 0; r < grid.dim; r++)
   {
      delete [] D1u[r];
      delete [] D1ucoef[r];
      delete [] D1uxcoef[r];
      delete [] D1uxxcoef[r];
   }
   delete [] D1u;
   delete [] D1ucoef;
   delete [] D1uxcoef;
   delete [] D1uxxcoef;
   free_matrix(D1jumpuxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);

   free_matrix(D2u,grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      for (s = 0; s < grid.dim; s++)
         free_matrix(D2ucoef[r][s],N,N,N);
      delete [] D2ucoef[r];
   }
   delete [] D2ucoef;
   free_matrix(D2uxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D2uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D2jumpuxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(jumpD1ucoef,N,N,N);
   free_matrix(jumpD1jumpuxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(jumpD2ucoef,N,N,N);

   for (r = 0; r < grid.dim; r++)
      free_matrix(uxcoef[r],N,N,N);
   delete [] uxcoef;
}

//
void cim345cond(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, 
            int *index, double***a, int gamma[][2], double ***S, PBData &pb, GridData &grid){

   double **G = matrix(2*grid.dim-1,2*grid.dim-1);

   char tempglobdist = globdist;

   globdist = 0;
   gmatrix(G, index, gamma, S, a, pb, grid);
   double estcond0 = condest(G, 6);


   globdist = 1;   
   gmatrix(G, index, gamma, S, a, pb, grid);
   double estcond1= condest(G, 6);

   globdist = (estcond0<estcond1)?0:1;
   
   free_matrix(G,2*grid.dim-1,2*grid.dim-1);

   if(a){
    cim345(A, b, Dusmall, buildsize, index, a, gamma, S, pb, grid);
   }else{
    cim345(A, b, Dusmall, buildsize, index, gamma, S, pb, grid);
   }
   
   globdist = tempglobdist;
  
}



//return a vector of all D2u scheme
char yescim5D2All(vector<double***> &D2ucoefvec, vector<double*> &D2uxcoefvec, vector<double*> &D2uxxcoefvec, vector<vector<int>> &offsetvec, int m, int n, 
               int *index, int mid, double ***S, GridData &grid)
{
   int i, r, s, thesign, N = 2*mid, ntheta = 8, signs[ntheta], sk[2], offset[2][2];
   int tindex[grid.dim], sindex[grid.dim];
   double theta0, theta, dtheta, themax, value;

   for (i = 0; i < grid.dim; i++)
      sindex[i] = mid;
   for (r = 0; r < grid.dim; r++)
      tindex[r] = index[r];

   thesign = 2*(evalarray(S,tindex) >= 0.0)-1;//1=outside, -1 indside
   dtheta = 2.0*M_PI/ntheta;
   theta0 = -3.0*M_PI/4.0;// start from bottom left corner
// computes signs of points surrounding node, start from (-1,-1) corner, couter-clockwise
   for (r = 0; r < ntheta; r++)
   {
      theta = theta0+r*dtheta;
      themax = max(fabs(cos(theta)),fabs(sin(theta)));
      offset[0][0] = round(cos(theta)/themax);
      offset[0][1] = round(sin(theta)/themax);
      tindex[m] = index[m]+offset[0][0];
      tindex[n] = index[n]+offset[0][1];
      signs[r] = 2*(evalarray(S,tindex) >= 0.0)-1;
   }
   tindex[m] = index[m];
   tindex[n] = index[n];

// looks for central differencing possibility around node
   for (r = 0; r < ntheta; r += 2)
      if ((thesign < 0)+(signs[r] < 0) != 1 && 
          (thesign < 0)+(signs[(r+2)%ntheta] < 0) != 1)
      {
         double ***D2ucoef = matrix(N,N,N);
         double *D2uxcoef = new double[3]();
         double *D2uxxcoef = new double[3]();

         theta = theta0+r*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[0][0] = round(cos(theta)/themax);
         offset[0][1] = round(sin(theta)/themax);
         theta = theta0+((r+2)%ntheta)*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[1][0] = round(cos(theta)/themax);
         offset[1][1] = round(sin(theta)/themax);
         sk[0] = offset[0][0]-offset[1][0];
         sk[1] = offset[0][1]-offset[1][1];
         if (sk[0] == 0)
         {
            sk[0] = offset[0][0];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[n] = -1.0/(sk[0]*grid.dx[m]);
         }
         else
         {
            sk[1] = offset[0][1];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[m] = -1.0/(sk[1]*grid.dx[n]);
         }
         sindex[m] = mid;
         sindex[n] = mid;

         offsetvec.push_back(vector<int>{offset[0][0],offset[0][1],offset[1][0],offset[1][1]});
         D2ucoefvec.push_back(D2ucoef);
         D2uxcoefvec.push_back(D2uxcoef);
         D2uxxcoefvec.push_back(D2uxxcoef);
         free_matrix(D2ucoef,N,N,N);
         delete [] D2uxcoef;
         delete [] D2uxxcoef;
      }

// looks for forward differencing possibility around node
   for (r = 0; r < ntheta; r++)
      if ((thesign < 0)+(signs[r] < 0) != 1 && 
          (thesign < 0)+(signs[(r+1)%ntheta] < 0) != 1)
      {
         double ***D2ucoef = matrix(N,N,N);
         double *D2uxcoef = new double[3]();
         double *D2uxxcoef = new double[3]();

         theta = theta0+r*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[0][0] = round(cos(theta)/themax);
         offset[0][1] = round(sin(theta)/themax);
         theta = theta0+((r+1)%ntheta)*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[1][0] = round(cos(theta)/themax);
         offset[1][1] = round(sin(theta)/themax);
         sk[0] = offset[0][0]-offset[1][0];
         sk[1] = offset[0][1]-offset[1][1];
         if (sk[0] == 0)
         {
            sk[0] = offset[0][0];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[n] = -1.0/(sk[0]*grid.dx[m]);
//            D2uxxcoef[n] = -0.5*sk[1]*grid.dx[n]/(sk[0]*grid.dx[m]);
            D2uxxcoef[n] = -0.5*(offset[0][1]+offset[1][1])*grid.dx[n]/
                            (sk[0]*grid.dx[m]);
         }
         else
         {
            sk[1] = offset[0][1];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[m] = -1.0/(sk[1]*grid.dx[n]);
//            D2uxxcoef[m] = -0.5*sk[0]*grid.dx[m]/(sk[1]*grid.dx[n]);
            D2uxxcoef[m] = -0.5*(offset[0][0]+offset[1][0])*grid.dx[m]/
                           (sk[1]*grid.dx[n]);
         }
         sindex[m] = mid;
         sindex[n] = mid;

         offsetvec.push_back(vector<int>{offset[0][0],offset[0][1],offset[1][0],offset[1][1]});
         D2ucoefvec.push_back(D2ucoef);
         D2uxcoefvec.push_back(D2uxcoef);
         D2uxxcoefvec.push_back(D2uxxcoef);
         free_matrix(D2ucoef,N,N,N);
         delete [] D2uxcoef;
         delete [] D2uxxcoef;
      }
   if (offsetvec.size() == 0){
      return 0;
   }
   return 1;
   
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

