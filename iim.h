#ifndef IIM_H
#define IIM_H

#include "grid.h"
#include "input.h"
#include "sparse.h"




struct ZLHeapElt
{
   int index[3];
   ZLHeapElt* child[2];
   ZLHeapElt* parent;
   int num;
   double val;
};

struct ZLHeapStruct
{
   ZLHeapElt *head;
   ZLHeapElt *tail;
   char ***mark;
   int dim;
};


void getiimjumps(double &up0, double &upum, double *uxp0, double **uxpuxm, 
                 double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm,
                 int *index, int rstar, int sk, double alpha, int thesign, 
                 double *normal, double ***S, PBData &pb, GridData &grid);
void getiimjumps(double &up0, double &upum, double *uxp0, double **uxpuxm, 
                 double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm,
                 double *x, int thesign, double ***S, PBData &pb, GridData &grid);
void getiimstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                     double up0, double upum, double *uxp0, double **uxpuxm, 
                     double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm, 
                     double ***S, PBData &pb, GridData &grid);
void getiimstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                     char yesC, double up0, double upum, double *uxp0, 
                     double **uxpuxm, double **uxxp0, double ***uxxpuxm, 
                     double ****uxxpuxxm, double ***S, PBData &pb, GridData &grid);
void getiimCstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                      double up0, double upum, double *uxp0, double **uxpuxm, 
                      double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm, 
                      double ***S, PBData &pb, GridData &grid);
void getiimgridstencilmx(double **A, int M, int N, double *x, int *gindex, int thesign, 
                         int **index, char yesC, double up0, double upum, double *uxp0, 
                         double **uxpuxm, double **uxxp0, double ***uxxpuxm, 
                         double ****uxxpuxxm, double ***S, PBData &pb, GridData &grid);
void iim(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S,
         PBData &pb, GridData &grid);
void iim(SparseElt2**** &A, double ***b, int *index, int add, double ***S, 
         char ***tube, PBData &pb, GridData &grid);
void iim(SparseElt2**** &A, double ***b, int *index, char yesC, int add, double ***S, 
         char ***tube, PBData &pb, GridData &grid);
void iimghost(SparseElt2**** &A, double ***b, int *index, char yesC, int add, 
              double ***S, char ***tube, PBData &pb, GridData &grid);
void iimC(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S,
          PBData &pb, GridData &grid);
void iimC(SparseElt2**** &A, double ***b, int *index, int add, double ***S, 
          char ***tube, PBData &pb, GridData &grid);
void addtoheap(ZLHeapStruct &heap, int *index, double val);
void fixheapeltup(ZLHeapStruct &heap, ZLHeapElt *fix);
ZLHeapElt* fixheapeltempty(ZLHeapStruct &heap, ZLHeapElt *fix);
void fixheapeltdelete(ZLHeapStruct &heap, ZLHeapElt *del);
ZLHeapElt* heapgoto(ZLHeapStruct &heap, int num);
void readZLHeap(ZLHeapStruct heap);
void getnearestgrid(int **index, int &N, double *x, int maxpts, double maxdist, 
                    char ***tube, GridData &grid);
double getdist(double *z, double *x, int thedim);
double bilinearinterp(double *x, double ***u, GridData &grid);
double weno6interpdirect(double *x, double ***u, GridData &grid);
void getallgrad(double ****grad, double ***S, GridData &grid);
int newtondir(double *x, double *x0, double *grad, double tol, double ***S, 
               double ****allgrad, GridData &grid);
int regulafalsidir(double *x, double *x0, double *grad, double tol, double ***S, 
                   GridData &grid);
void getnearestinterface(double *x, int *index, double ***S, GridData &grid) ;


void testZL(char yesC, int add, double ***S, char ***tube, PBData &pb, 
            GridData &grid);
void testZLatx(double *theerr, double *x, double thesign, char yesC, int add, 
               double ***S, char ***tube, PBData &pb, GridData &grid);

// used in IIM
void getinterfaceinfo(double *normal, double *tangent1, double *tangent2, double **Dn,
                      double &tau, double *Dtau, double **D2tau, double &sigma, 
                      double *Dsigma, double &jumpfe, double *x, double ***S,
                      PBData &pb, GridData &grid);


#endif