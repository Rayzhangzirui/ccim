#include "global.h"

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
int SURFOPT;


// #define FIXBANANA // surgical fix for banana shape at
bool globwritemx = false;//write coupling matrix for analysis
bool globwriteerr = false;//write error in u and Du at each grid point

