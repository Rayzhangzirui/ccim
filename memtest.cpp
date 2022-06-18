#include "tryvn.h"
#include "helper.h"

using namespace std;

extern double globtime;
extern int runStep;

extern int SURFOPT;
extern char globbiasstat;
extern int globtesta;
extern int globsmall;
extern int GRIDNUM;
extern int globtestnum;
extern int globperturb;
extern int globcim;

extern int globlinsolve;
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
extern bool globwriteerr;

extern ofstream outfile_mx;
extern ofstream outfile_uerr;
extern ofstream outfile_Duerr;

extern double RADIUS;

int main(int argc, char* argv[])
{

	cmdLine(argc, argv);
	
	SparseAMGMx P;
	P.row = new SparseAMGElt*[21466890];
	return 0;
}




