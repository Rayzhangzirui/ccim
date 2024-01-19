#include "global.h"
#include <cstdlib>
#include <stdio.h>
#include <getopt.h>
#include <iostream>

using namespace std;


double globerr = 0.0;
double globerrvec[3], globerrvec2[3], globerrvec3[3], globerrvec4[3]; 
int globerrmx[3][3]; // glboerrmx[r] index of max err u_rr

char globorder = 2, globintorder = 3;

// if globdist = 1, use getcim345D2udist, else getcim345D2u
// note globdist == 1 means cim 3 po ints may use cim 5 differencing
char globdist = 0, globdistvar = 1;
char globdirectD2 = 1; //used in getcim345D2u, use cross derivative at nbr po int in the same side
char globbiasstat = 1;// when approximating mixed derivative, prefer larger status?
char globcheck = 1;// check Du DDu
char globnorm = 2; // norm for residual
char globheap = 1; // used in AMG
char globoldwrongcim4 = 0, globoldwrongcim5 = 0,
globoldcim5 = 0, globdebug = 0;


char globheat = 1;//getheatvn, usual way to calculate vn
char globexactmotion = 0, globrk = 0, globcim1order = 1, globdtdx2 = 1;

char globexactvn = 0; //use exact vn
char globregfalweno = 1;

int globGSsmooth = 0;
int globheatsmooth = 0;
double globheatsmoothtime = 0.0;

int globtestnum = 0;

char globallcim345 = 1; //use cim345 for all 2nd order point
int globcim = 345;
// 1: cim 1 
// 2: cim 2 and cim 1
// 4: cim 3 and cim 4
// 5: cim 3 and cim 5
// 345: cim 3 and combination of cim 4 and cim 5
// 6: cim based on estimated condition number 
// 0: iim

double tollinsolve = 1e-8; // tolerance of linear solver

int globlinsolve = 1;
// 0: BICGSTAB
// 1: BICGSTAB with ILU preconditioning
// 2: AMG
// 3: GMRES
// 4: hypre amg
// 5: hypre bicgstab
int globperturb = 0;
// 0: CIM perturbation
// n > 0: pertrubation using status for maximum of n steps 
// -1: no perturbation

int globsmall = 1;
// 0: not storing info for Du
// 1: storing info for Du in local index
// 2: storing info for Du in global index

int GRIDNUM = 30;//global grid number
int eindex[3] = {-1,-1,-1};
double EPSILONP = 80.0;
double EPSILONM = 2.0;
int globtesta = 0; // a term
int SURFOPT = 0;


// #define FIXBANANA // surgical fix for banana shape at
bool globwritemx = false;//write coupling matrix for analysis
bool globwriteerr = false;//write error in u and Du at each grid point


// motion
double RADIUS = 0.5;
int runStep = 2000;
double globtime = 0.1;

// protein surface
double ETA = 30.0;
// files to store protein surface
string PROTEINDIR = "/Users/zziruiadmin/projects/ccimproj/ccim/protein_surf/"; 
string PROJECTDIR = "/Users/zziruiadmin/projects/ccimproj/ccim/";
string PROTEINFILE = "mdm2.txt";


bool globgridperturb = false;

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
            EPSILONM = atof(optarg);
         break;

         case 'P':
            EPSILONP = atof(optarg);
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
               "-s --SURFOPT: 0 = circle, 1 = wierd torus, 11 = ellipsoid, 12 = donut, 13 = banana, 14 = 8 balls, 15 = peanut, 16 = popcorn, 20 = protien surface\n"
               "-a --globtesta, 0 = without a term, 1 = with a term\n"
               "-l --globlinsolve, 0 = BICGSTAB, 1 = ILU BICGSTAB, 2 = AMG, 4 = hypre-amg, 5 = hypre-ilu-bicgsta\n"
               "-k --globcheck, 1 = check convergence, 0 = don't check convergence, used for motion computation\n"
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
               "-R --radius, radius of circle\n"
               "-T --globgridperturb, 1 =  randomly shift grid\n"
               "-o --tollinsolve, tolerance of linear system solver\n"
               "-f --PROTEINFILE\n");
               
             exit(-1);
      }
   }
   for (; optind < argc; optind++)
      printf ("Non-option argument %s\n", argv[optind]);
   
   if (globtestnum == 4){
      printf("test with a\n");
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