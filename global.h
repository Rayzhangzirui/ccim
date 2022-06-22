#ifndef GLOBAL_H
#define GLOBAL_H
#include <string>
#include <fstream>


extern double globerr;
extern double globerrvec[3], globerrvec2[3], globerrvec3[3], globerrvec4[3]; 
extern int globerrmx[3][3]; // glboerrmx[r] index of max err u_rr

extern char globorder, globintorder;
// if globdist = 1, use getcim345D2udist, else getcim345D2u
// note globdist == 1 means cim 3 poextern ints may use cim 5 differencing
extern char globdist, globdistvar;
extern char globdirectD2 ; //used in getcim345D2u, use cross derivative at nbr point in the same side
extern char globbiasstat ;// when approximating mixed derivative, prefer larger status?
extern char globcheck ;// check Du DDu
extern char globnorm ; // norm for residual
extern char globheap ; // used in AMG
extern char globoldwrongcim4 , globoldwrongcim5,
     globoldcim5, globdebug;


extern char globheat;//getheatvn, usual way to calculate vn
extern char globexactmotion, globrk, globcim1order, globdtdx2;

extern char globexactvn; //use exact vn
extern char globregfalweno;

extern int globGSsmooth;
extern int globheatsmooth;
extern double globheatsmoothtime;

extern int globtestnum;

extern char globallcim345; //use cim345 for all 2nd order poextern int
extern int globcim;
// 1: cim 1 
// 2: cim 2 and cim 1
// 4: cim 3 and cim 4
// 5: cim 3 and cim 5
// 345: cim 3 and combination of cim 4 and cim 5
// 0: iim

extern double tollinsolve; // tolerance of linear solver

extern int globlinsolve;
// 0: BICGSTAB
// 1: BICGSTAB with ILU preconditioning
// 2: AMG
// 3: GMRES
// 4: hypre amg
// 5: hypre bicgstab
extern int globperturb;
// 0: CIM perturbation
// n > 0: pertrubation using status for maximum of n steps 
// -1: no perturbation
extern int globsmall;
// 0: not storing info for Du
// 1: storing info for Du in local index
// 2: storing info for Du in global index
extern int GRIDNUM;//global grid number
extern int eindex[3];
extern double EPSILONP;
extern double EPSILONM;
extern int globtesta; // a term
extern int SURFOPT;


// #define FIXBANANA // surgical fix for banana shape at
extern bool globwritemx;//write coupling matrix for analysis
extern bool globwriteerr;//write error in u and Du at each grid point

// motion
extern double RADIUS;
extern int runStep;
extern double globtime;

// protein surface
extern double ETA ;
extern std::string PROTEINDIR ;
extern std::string PROJECTDIR ;
extern std::string PROTEINFILE ;


extern bool globgridperturb;



extern std::ofstream outfile_mx;
extern std::ofstream outfile_uerr;
extern std::ofstream outfile_Duerr;


void cmdLine(int argc, char *argv[]);

#endif