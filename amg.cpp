#include "amg.h"
#include "numerics.h"
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

extern int globnorm;
extern int globheap;

// argument theta = 1
void startamg(int *strength, SparseAMGMx &A, double theta)
{
   int i, j;
   double themax, tol = 1.0-1.0e-14;

   for (i = 0; i < A.n; i++)
      strength[i] = 0;

   for (i = 0; i < A.n; i++)
   {
      themax = 0.0;
      for (j = 0; j < A.nc[i]; j++)
         if ((A.row[i][j]).roc != i && themax < -(A.row[i][j]).val[0]) 
            themax = -(A.row[i][j]).val[0];
      for (j = 0; j < A.nc[i]; j++)
         if (-(A.row[i][j]).val[0] > tol*theta*themax) 
         {
            (A.row[i][j]).strong[0] = 1;
            (strength[(A.row[i][j]).roc])++;
         }
         else
            (A.row[i][j]).strong[0] = 0;
   }
}

void startamg(int *strength, int &strhead, int* &strnext, int* &strprev, 
              SparseAMGMx &A, double theta)
{
   int i, j;
   double themax, tol = 1.0-1.0e-14;

   strhead = 0;
   strnext = new int[A.n]; 
   strprev = new int[A.n];
   for (i = 0; i < A.n; i++)
   {
      strength[i] = 0;
      strnext[i] = i+1;
      strprev[i] = i-1;
   }

   for (i = 0; i < A.n; i++)
   {
      themax = 0.0;
      for (j = 0; j < A.nc[i]; j++)
         if ((A.row[i][j]).roc != i && themax < -(A.row[i][j]).val[0])
            themax = -(A.row[i][j]).val[0];
      for (j = 0; j < A.nc[i]; j++)
         if (-(A.row[i][j]).val[0] > tol*theta*themax)
         {
            (A.row[i][j]).strong[0] = 1;
            (strength[(A.row[i][j]).roc])++;
         }
         else
            (A.row[i][j]).strong[0] = 0;
   }
}

void startamg(int *strength, StrHeapStruct &heap, SparseAMGMx &A, double theta)
{
   int i, j;
   double themax, tol = 1.0-1.0e-14;

   for (i = 0; i < A.n; i++)
      strength[i] = 0;

   for (i = 0; i < A.n; i++)
   {
      themax = 0.0;
      for (j = 0; j < A.nc[i]; j++)
         if ((A.row[i][j]).roc != i && themax < -(A.row[i][j]).val[0]) // get max of negative of non-diagonal
            themax = -(A.row[i][j]).val[0];
      for (j = 0; j < A.nc[i]; j++)
         if (-(A.row[i][j]).val[0] > tol*theta*themax)
         {
            (A.row[i][j]).strong[0] = 1; // variable A.row[i][j].roc strongly depends on variable i 
            (strength[(A.row[i][j]).roc])++; 
         }
         else
            (A.row[i][j]).strong[0] = 0;
   }

   heap.value = strength;
   heap.index = new int[A.n];
   heap.loc = new int[A.n];
   heap.num = 0;
   for (i = 0; i < A.n; i++)
      addtoheap(heap,i);
}

void amgfirstpass(int *coarse, int *strength, SparseAMGMx &A)
{
   int i, j, k, l, maxval, maxloc;

   for (i = 0; i < A.n; i++)
      coarse[i] = strength[i];

   maxval = -1; 
   maxloc = -1;
   for (i = 0; i < A.n; i++)
      if (coarse[i] > maxval)
      {
         maxval = coarse[i];
         maxloc = i;
      }

   while (maxval >= 0)
   {
      coarse[maxloc] = -2;
      for (i = 0; i < A.nr[maxloc]; i++)
      {
         k = (A.col[maxloc][i]).roc;
         if ((A.col[maxloc][i]).strong[0] && coarse[k] >= 0) 
         {
            coarse[k] = -1;
            for (j = 0; j < A.nc[k]; j++)
            {
               l = (A.row[k][j]).roc;
               if ((A.row[k][j]).strong[0] && coarse[l] >= 0)
                  (coarse[l])++;
            }
         }
      }

      maxval = -1;
      maxloc = -1;
      for (i = 0; i < A.n; i++)
         if (coarse[i] > maxval)
         {
            maxval = coarse[i];
            maxloc = i;
         }
   }

   for (i = 0; i < A.n; i++)
      coarse[i] = -(coarse[i]+1);
}

void amgfirstpass(int *coarse, int *strength, int &strhead, int* &strnext, 
                  int* &strprev, SparseAMGMx &A)
{
   int i, j, k, l, maxval, maxloc;

   for (i = 0; i < A.n; i++)
      coarse[i] = strength[i];

   while (strhead < A.n)
   {
      maxval = -1;
      maxloc = -1;
      for (i = strhead; i < A.n; i = strnext[i])
         if (coarse[i] > maxval)
         {
            maxval = coarse[i];
            maxloc = i;
         }
      coarse[maxloc] = -2;
      if (strprev[maxloc] >= 0)
         strnext[strprev[maxloc]] = strnext[maxloc];
      else
         strhead = strnext[maxloc];
      if (strnext[maxloc] < A.n)
         strprev[strnext[maxloc]] = strprev[maxloc];

      for (i = 0; i < A.nr[maxloc]; i++)
      {
         k = (A.col[maxloc][i]).roc;
         if ((A.col[maxloc][i]).strong[0] && coarse[k] >= 0)
         {
            coarse[k] = -1;
            if (strprev[k] >= 0)
               strnext[strprev[k]] = strnext[k];
            else
               strhead = strnext[k];
            if (strnext[k] < A.n)
               strprev[strnext[k]] = strprev[k];

            for (j = 0; j < A.nc[k]; j++)
            {
               l = (A.row[k][j]).roc;
               if ((A.row[k][j]).strong[0] && coarse[l] >= 0)
                  (coarse[l])++;
            }
         }
      }
   }

   for (i = 0; i < A.n; i++)
      coarse[i] = -(coarse[i]+1);

   delete [] strnext;
   delete [] strprev;
}

void amgfirstpass(int *coarse, int *strength, StrHeapStruct &heap, SparseAMGMx &A)
{
   int i, j, k, l, maxval, maxloc;

   for (i = 0; i < A.n; i++)
      coarse[i] = strength[i];

   while (heap.num > 0)
   {
      maxloc = heap.index[0];
      coarse[maxloc] = -2; // mark as C-point
      popfromheap(heap,0);

      for (i = 0; i < A.nr[maxloc]; i++)
      {
         k = (A.col[maxloc][i]).roc;
         if ((A.col[maxloc][i]).strong[0] && coarse[k] >= 0) // if variable (A.col[maxloc][i]).roc strongly depends on i and not yet visitied
         {
            coarse[k] = -1; // mark is F-point
            popfromheap(heap,heap.loc[k]);

            for (j = 0; j < A.nc[k]; j++) //For each new F-point, the weights of its neighbors are increased by 1
            {
               l = (A.row[k][j]).roc;
               if ((A.row[k][j]).strong[0] && coarse[l] >= 0)
               {
                  (coarse[l])++;
                  fixheap(heap,heap.loc[l]);
               }
            }
         }
      }
   }

   for (i = 0; i < A.n; i++)
      coarse[i] = -(coarse[i]+1); // -1 -> 0, -2 -> 1

   delete [] heap.index; // do not delete heap.value = strength
   delete [] heap.loc;
}
// An Introduction to Algebraic Multigrid, R. D. Falgout
// The original C-AMG interpolation scheme (described
// below) requires each pair of stronglyconnected
// F-points to be strongly connected to a
// common C-point. The second pass of the coarsening
// algorithm searches for F-point pairs that do not
// satisfy this requirement, and changes one of them to
// a C-point
void amgsecondpass(int *coarse, SparseAMGMx &A)
{
   int i, j, k, l, r, s, count1, count2;
   char yescriteria;

   for (i = 0; i < A.n; i++)
   {
      if (!(coarse[i])) 
         for (j = 0; j < A.nc[i]; j++) 
         {
            k = (A.row[i][j]).roc;
            if ((A.row[i][j]).strong[0] && !(coarse[k])) 
            {
// points i and k both fine and i influenced by k
               yescriteria = 0;
               for (l = 0; l < A.nc[k] && !yescriteria; l++)
               {
                  r = (A.row[k][l]).roc; // r is connected with k
                  if ((A.row[k][l]).strong[0] && coarse[r])
                  {
// r is coarse and k influenced by r
// check i is also influenced by r
                     for (s = 0; s < A.nc[i] && 
                                 (!((A.row[i][s]).strong[0]) || (A.row[i][s]).roc != r); 
                          s++);
                     if (s < A.nc[i])
                        yescriteria = 1; // found common C-point
                  }
               }
               if (!yescriteria)
               {
                  count1 = 0;
                  for (l = 0; l < A.nc[i]; l++)
                     if ((A.row[i][l]).strong[0] && coarse[(A.row[i][l]).roc])
                        count1++;
                  for (l = 0; l < A.nr[i]; l++)
                     if ((A.col[i][l]).strong[0] && coarse[(A.col[i][l]).roc])
                        count1++;
                  count2 = 0;
                  for (l = 0; l < A.nc[k]; l++)
                     if ((A.row[k][l]).strong[0] && coarse[(A.row[k][l]).roc])
                        count2++;
                  for (l = 0; l < A.nr[k]; l++)
                     if ((A.col[k][l]).strong[0] && coarse[(A.col[k][l]).roc])
                        count2++;

                  if (count1 <= count2)
                     coarse[i] = 1;
                  else
                     coarse[k] = 1;
               }
            }
         }
   }
}
// coarse and strength are the same array
void makecoarsemx(int *coarse, int *strength, SparseAMGMx &A)
{
   cout << "   first pass" << endl;
   amgfirstpass(coarse,strength,A);
   cout << "   second pass" << endl;
   amgsecondpass(coarse,A);

/*
   for (int i = 0; i < A.n; i++)
      for (int r = 0; r < A.nc[i]; r++)
         cout << (A.row[i][r]).roc << " " << (A.row[i][r]).val[0] << " "
              << (int) (A.row[i][r]).strong[0] << endl;
   for (int j = 0; j < A.m; j++)
      for (int r = 0; r < A.nr[j]; r++)
         cout << (A.col[j][r]).roc << " " << (A.col[j][r]).val[0] << " "
              << (int) (A.col[j][r]).strong[0] << endl;
*/
}

void makecoarsemx(int *coarse, int *strength, int &strhead, int *strnext, int *strprev, 
                  SparseAMGMx &A)
{
   clock_t cstart1, cend1, cstart2, cend2;

   cout << "   first pass" << endl;
   cstart1 = clock ();
   amgfirstpass(coarse,strength,strhead,strnext,strprev,A);
   cend1 = clock ();
   cstart2 = clock ();
   cout << "   second pass" << endl;
   amgsecondpass(coarse,A);
   cend2 = clock ();
   cout << "time = " << (double) (cend1-cstart1)/CLOCKS_PER_SEC << " and " 
        << (double) (cend2-cstart2)/CLOCKS_PER_SEC << endl;
}

void makecoarsemx(int *coarse, int *strength, StrHeapStruct &heap, SparseAMGMx &A)
{
//   clock_t cstart1, cend1, cstart2, cend2;

//   cout << "   first pass" << endl;
//   cstart1 = clock ();
   amgfirstpass(coarse,strength,heap,A);
//   cend1 = clock ();
//   cstart2 = clock ();
//   cout << "   second pass" << endl;
   amgsecondpass(coarse,A);
//   cend2 = clock ();
//   cout << "time = " << (double) (cend1-cstart1)/CLOCKS_PER_SEC << " and " 
//        << (double) (cend2-cstart2)/CLOCKS_PER_SEC << endl;
}

char yesinfluencedby(int i, int j, SparseAMGMx &A)
{
   int k;

   for (k = 0; k < A.nc[i]; k++)
      if ((A.row[i][k]).strong[0] && (A.row[i][k]).roc == j)
         return 1;

   return 0;
}



SparseAMGMx createP(int *coarse, SparseAMGMx &A)
{
   int i, j, k, l, r, s, count, order[A.n];
   double denom, tempdenom, *temp;
   SparseElt **Pt, *current, *current2;
   SparseAMGMx P;
 
   P.n = A.n;
   P.row = new SparseAMGElt*[P.n];
   P.nc = new int[P.n];
   count = 0; // total number of C-point, dimension of P is n x count
   for (i = 0; i < A.n; i++)
      if (coarse[i])
      {
         order[i] = count;
         count++;
      }
   P.m = count;
   P.col = new SparseAMGElt*[P.m];
   P.nr = new int[P.m];

   Pt = new SparseElt*[P.m];
   for (i = 0; i < P.m; i++)
      Pt[i] = NULL;

   for (i = 0; i < P.n; i++)
      if (coarse[i])
      {
         P.nc[i] = 1;
         P.row[i] = new SparseAMGElt[P.nc[i]];
         (P.row[i][0]).roc = order[i];
         (P.row[i][0]).val = new double[1];
         (P.row[i][0]).val[0] = 1.0;
         (P.row[i][0]).strong = NULL;
         sparseorder((P.row[i][0]).roc, i, Pt, (P.row[i][0]).val[0]);
      }
      else
      {
// i is fine
         count = 0;
         denom = 0.0;
         for (j = 0; j < A.nc[i]; j++)
            if ((A.row[i][j]).strong[0] && coarse[(A.row[i][j]).roc]) //count strongly connected C-point
               count++;
            else if (!((A.row[i][j]).strong[0]))
// assuming diagonal entry is a weak connection
               denom += (A.row[i][j]).val[0];
         P.nc[i] = count;
         P.row[i] = new SparseAMGElt[P.nc[i]];
         temp = new double[P.nc[i]]; //coefficient for each strongly connnected C-points

         count = 0;
         for (j = 0; j < A.nc[i]; j++)
         {
            k = (A.row[i][j]).roc;
// k is coarse and i is influenced by k
            if ((A.row[i][j]).strong[0] && coarse[k])
            {
               (P.row[i][count]).roc = order[k];
               (P.row[i][count]).val = new double[1];
               (P.row[i][count]).val[0] = -(A.row[i][j]).val[0];
               (P.row[i][count]).strong = NULL;
               count++;
            }
         }
         for (j = 0; j < A.nc[i]; j++)
         {
            k = (A.row[i][j]).roc;
// k is fine and i is influenced by k
            if ((A.row[i][j]).strong[0] && !(coarse[k]))
            {
               tempdenom = 0.0;
               for (l = 0; l < P.nc[i]; l++)
                  temp[l] = 0.0;
               for (l = 0; l < A.nc[k]; l++)
               {
                  r = (A.row[k][l]).roc;
// r is coarse and k is influenced by r
                  if ((A.row[k][l]).strong[0] && coarse[r])
                  {
                     count = 0;
                     for (s = 0; s < A.nc[i] && count >= 0; s++)
                        if ((A.row[i][s]).strong[0] && coarse[(A.row[i][s]).roc]) // (A.row[i][s]) is coarse and influence i
                        {
                           if ((A.row[i][s]).roc != r)
                              count++;
                           else
                           {
// i is also influenced by r
// add up contributions of k in terms of r
                              temp[count] += (A.row[i][j]).val[0]*
                                             ((A.row[k][l]).val[0]); 
                              tempdenom += (A.row[k][l]).val[0];
                              count = -1;//break s-loop
                           }
                        }
                  }
               }
               for (l = 0; l < P.nc[i]; l++)
                  (P.row[i][l]).val[0] -= temp[l]/tempdenom;
            }
         }
         delete[] temp;
         for (j = 0; j < P.nc[i]; j++)
         {
            (P.row[i][j]).val[0] /= denom;
//            if (Pt[0] != NULL)
//               cout << (*(Pt[0])).cindex << " " << i << " " << j << " " 
//                    << (P.row[i][j]).roc << endl;
            sparseorder((P.row[i][j]).roc,i,Pt,(P.row[i][j]).val[0]);
//            if (Pt[0] != NULL)
//               cout << (*(Pt[0])).cindex << " " << i << " " << j << " " 
//                    << (P.row[i][j]).roc << endl;
         }
      }

//   for (i = 0; i < P.n; i++)
//      for (j = 0; j < P.nc[i]; j++)
//         cout << i << " " << (P.row[i][j]).roc << " " << (P.row[i][j]).val[0] << endl;
//   cout << P.m << endl;
//   for (j = 0; j < P.m; j++)
//      for (current = Pt[j]; current != NULL; current = (*current).next)
//      {
//         cout << (*current).val << endl;
//         cout << (*current).cindex << endl;
//         cout << (*current).cindex << " " << j << " " << (*current).val << endl;
//      }

   for (j = 0; j < P.m; j++)
   {
      count = 0;
      for (current = Pt[j]; current != NULL; current = (*current).next)
         count++;
      P.nr[j] = count;
      P.col[j] = new SparseAMGElt[P.nr[j]];
      count = 0;
      for (current = Pt[j]; current != NULL; 
           current2 = current, current = (*current).next, delete current2)
      {
         (P.col[j][count]).roc = (*current).cindex;
         i = (P.col[j][count]).roc;
         for (k = 0; k < P.nc[i] && (P.row[i][k]).roc != j; k++);
         (P.col[j][count]).val = (P.row[i][k]).val;
         (P.col[j][count]).strong = NULL;
         count++;
      }
   }

   delete[] Pt;

   return P;
}

double getdotprod(SparseAMGElt *v, int numv, SparseAMGElt *w, int numw)
{
// assumes ordered
   int i, j;
   double value = 0.0;

   for (i = 0,j = 0; i < numv && j < numw; )
      if ((v[i]).roc == (w[j]).roc)
      {
         value += (v[i]).val[0]*((w[j]).val[0]);
         i++;  
         j++;
      }
      else if ((v[i]).roc < (w[j]).roc)
         i++;
      else
         j++;

   return value;
}

double getdotprod(SparseAMGElt *v, int numv, SparseElt *w)
{
// assumes ordered
   int i;
   double value = 0.0;
   SparseElt *current;

   for (i = 0, current = w; i < numv && current != NULL; )
      if ((v[i]).roc == (*current).cindex)
      {
         value += (v[i]).val[0]*((*current).cindex);
         i++;  
         current = (*current).next;
      }
      else if ((v[i]).roc < (*current).cindex)
         i++;
      else
         current = (*current).next;

   return value;
}

double getdotprod(SparseAMGElt *v, int numv, double *w)
{
// assumes ordered
   int i;
   double value = 0.0;

   for (i = 0; i < numv; i++)
      value += (v[i]).val[0]*w[(v[i]).roc];

   return value;
}

SparseAMGMx multPtP(SparseAMGMx &A, SparseAMGMx &P)
{
   int i, j, k, count;
   SparseAMGMx B;
   SparseElt *temp = NULL, *temp2 = NULL, *current, **C;
   double value, tol = 1.0e-14;

   C = new SparseElt*[P.m];
   for (i = 0; i < P.m; i++)
      C[i] = NULL;

   B.n = P.m;
   B.m = P.m;
   B.row = new SparseAMGElt*[B.n];
   B.nc = new int[B.n];
   B.col = new SparseAMGElt*[B.m];
   B.nr = new int[B.m];

   for (j = P.m-1; j >= 0; j--)
   {
      for (k = A.n-1; k >= 0; k--)
      {
         value = getdotprod(A.row[k],A.nc[k],P.col[j],P.nr[j]);
         if (fabs(value) > tol)
         {
            current = temp;
            temp = new SparseElt;
            (*temp).val = value;
            (*temp).cindex = k;
            (*temp).next = current;
         }
      }
      count = 0;
      for (i = P.m-1; i >= 0; i--)
      {
         value = getdotprod(P.col[i],P.nr[i],temp);
         if (fabs(value) > tol)
         {
            current = temp2;
            temp2 = new SparseElt;
            (*temp2).val = value;
            (*temp2).cindex = i;
            (*temp2).next = current;
            count++;
         }
         
      }
      for (current = temp; current != NULL; 
           temp = current, current = (*current).next, delete temp);
      temp = NULL;
      B.col[j] = new SparseAMGElt[count];
      count = 0;
      for (current = temp2; current != NULL; 
           temp2 = current, current = (*current).next, delete temp2)
      {
         (B.col[j][count]).roc = (*current).cindex;
         (B.col[j][count]).val = new double[1];
         (B.col[j][count]).val[0] = (*current).val;
         (B.col[j][count]).strong = new char[1];
         sparseorder((B.col[j][count]).roc,j,C,(B.col[j][count]).val[0]);
         count++;
      }
      temp2 = NULL;
   }

   for (i = 0; i < B.n; i++)
   {
      count = 0;
      for (current = C[i]; current != NULL; current = (*current).next)
         count++;
      B.row[i] = new SparseAMGElt[count];
      count = 0;
      for (current = C[i]; current != NULL; 
           temp = current, current = (*current).next, delete temp)
      {
         (B.row[i][count]).roc = (*current).cindex;
         j = (B.row[i][count]).roc;
         for (k = 0; k < B.nr[j] && (B.col[j][k]).roc != i; k++);
         (B.row[i][count]).val = (B.col[j][k]).val;
         (B.row[i][count]).strong = (B.col[j][k]).strong;
         count++;
      }
   }
   delete [] C;

   return B;
}

SparseAMGMx multPtP2(SparseAMGMx &A, SparseAMGMx &P)
{
   int i, j, k, r, s, count;
   SparseAMGMx B;
   SparseElt *temp = NULL, *temp2 = NULL, *current, **C;
   double value, tol = 1.0e-14;

   C = new SparseElt*[P.m];
   for (i = 0; i < P.m; i++)
      C[i] = NULL;

   B.n = P.m;
   B.m = P.m;
   B.row = new SparseAMGElt*[B.n];
   B.nc = new int[B.n];
   B.col = new SparseAMGElt*[B.m];
   B.nr = new int[B.m];
   
   for (j = P.m-1; j >= 0; j--)
   {
// jth column of P
      temp = NULL;
      for (r = 0; r < P.nr[j]; r++)
      {
         k = (P.col[j][r]).roc;
// use elements, such as kth, in jth column as weights in linear combination
         for (s = 0; s < A.nr[k]; s++)
// linear combination of columns of A, such as kth column
            sparseorder((A.col[k][s]).roc,temp,
                        (P.col[j][r]).val[0]*((A.col[k][s]).val[0]));
      }
// use elements in column for linear combination
      temp2 = NULL;
      for (current = temp; current != NULL; 
           temp = current, current = (*current).next, delete temp)
         if (fabs((*current).val) > tol)
         {
// use kth element
            k = (*current).cindex;
            for (s = 0; s < P.nc[k]; s++)
// linear combination of P', such as kth column
               sparseorder((P.row[k][s]).roc,temp2,
                           (*current).val*((P.row[k][s]).val[0]));
         }

      count = 0;
      for (current = temp2; current != NULL; current = (*current).next)
         if (fabs((*current).val) > tol)
            count++;
      B.nr[j] = count;
      B.col[j] = new SparseAMGElt[B.nr[j]];
      count = 0;
      for (current = temp2; current != NULL; 
           temp2 = current, current = (*current).next, delete temp2)
         if (fabs((*current).val) > tol)
         {
            (B.col[j][count]).roc = (*current).cindex;
            (B.col[j][count]).val = new double[1];
            (B.col[j][count]).val[0] = (*current).val;
            (B.col[j][count]).strong = new char[1];
            sparseorder((B.col[j][count]).roc,j,C,(B.col[j][count]).val[0]);
            count++;
         }
   }

   for (i = 0; i < B.n; i++)
   {
      count = 0;
      for (current = C[i]; current != NULL; current = (*current).next)
         count++;
      B.nc[i] = count;
      B.row[i] = new SparseAMGElt[B.nc[i]];
      count = 0;
      for (current = C[i]; current != NULL; 
           temp = current, current = (*current).next, delete temp)
      {
         (B.row[i][count]).roc = (*current).cindex;
         j = (B.row[i][count]).roc;
         for (r = 0; r < B.nr[j] && (B.col[j][r]).roc != i; r++);
         (B.row[i][count]).val = (B.col[j][r]).val;
         (B.row[i][count]).strong = (B.col[j][r]).strong;
         count++;
      }
   }
   delete [] C;

   return B;
}

void clearsparse(SparseAMGMx &A)
{
   int i, j;
   
   for (i = 0; i < A.n; i++)
   {
      for (j = 0; j < A.nc[i]; j++)
      {
         delete (A.row[i][j]).val;
         delete (A.row[i][j]).strong;
      }
      delete [] A.row[i];
   }
   if (A.n > 0)
   {
      delete [] A.row;
      delete [] A.nc;
   }

   for (j = 0; j < A.m; j++)
      delete [] A.col[j];
   if (A.m > 0)
   {
      delete [] A.col;
      delete [] A.nr;
   }

   A.n = 0;
   A.m = 0;
}

void mxvecmult(double *a, SparseAMGMx A, double *b)
{
   int i;

   for (i = 0; i < A.n; i++)
      a[i] = getdotprod(A.row[i],A.nc[i],b);
}

void mxtvecmult(double *a, SparseAMGMx A, double *b)
{
   int i;

   for (i = 0; i < A.m; i++)
      a[i] = getdotprod(A.col[i],A.nr[i],b);
}

void interpolate(SparseAMGMx &B, double* &c, SparseAMGMx &P, int *coarse, 
                 SparseAMGMx &A, double *b)
{
   int i;

   cout << "   createP" << endl;
   P = createP(coarse,A);
//   checkmx(P);
   cout << "   multPtP" << endl;
   B = multPtP2(A,P);
//   checkmx(B);
   c = new double[P.m];
   for (i = 0; i < P.m; i++)
      c[i] = getdotprod(P.col[i],P.nr[i],b);
}

void interpolate(SparseAMGMx &B, SparseAMGMx &P, int *coarse, SparseAMGMx &A)
{
   int i;

   P = createP(coarse,A);
   B = multPtP2(A,P);
}

void gaussseidel(double *x, SparseAMGMx &A, double *b, double *x0, int nu)
{
   int i, j, step;
   double diag;

   for (i = 0; i < A.n; i++)
      x[i] = x0[i];

   for (step = 1; step <= nu; step++)
      for (i = 0; i < A.n; i++)
      {
         x[i] = b[i];
         for (j = 0; j < A.nc[i]; j++)
            if ((A.row[i][j]).roc != i)
               x[i] -= (A.row[i][j]).val[0]*x[(A.row[i][j]).roc];
            else
               diag = (A.row[i][j]).val[0];
         x[i] /= diag;
      }
}

double Vcycle(double *x, SparseAMGMx &A, double *b, double *x0, int nu1, int nu2, 
              int small)
{
   int i, j, *coarse = new int[A.n];
   SparseAMGMx Acoarse, P;
   double rnorm, theta = 1.0, *rcoarse = NULL, *xcoarse = NULL, *r = new double[A.n];
   
/*
   delta2d.precision(16);

   delta2d << A.n << " " << A.m << endl;
   delta2d << "b" << endl;
   for (i = 0; i < A.n; i++)
      delta2d << scientific << b[i] << " ";
   delta2d << endl;
   delta2d << "A" << endl;
   outputAMGmx(delta2d,A);
*/
//   checkmx(A);
   cout << A.n << " " << A.m << endl;
//   peekmx(A);
//   getchar();

   cout << "in vcycle" << endl;
   for (i = 0; i < A.n; i++)
      x[i] = x0[i];

   if (A.n <= small) 
   {
      cout << "small enough " << A.n << endl;
      double **B = matrix(A.n-1,A.m-1); 
      int *PLR = new int[A.n], *PLC = new int[A.m];
      
      for (i = 0; i < A.n; i++) 
         for (j = 0; j < A.m; j++) 
            B[i][j] = 0.0;
      for (i = 0; i < A.n; i++)
         for (j = 0; j < A.nc[i]; j++) 
            B[i][(A.row[i][j]).roc] = (A.row[i][j]).val[0];
      gecp0(B,PLR,PLC,B,A.n-1,A.m-1);
      forwardbacksub0(x,b,B,PLR,PLC,A.n-1);

      free_matrix(B,A.n-1,A.m-1);
      delete [] PLR;
      delete [] PLC;
   }
   else 
   {
      cout << "gauss seidel" << endl;
      gaussseidel(x,A,b,x,nu1);

// calculate residual
      cout << "residual" << endl;
      mxvecmult(r,A,x);
      for (i = 0; i < A.n; i++)
         r[i] = b[i]-r[i];

      cout << "here" << endl;
//      int strhead, strnext[A.n], strprev[A.n];
// make coarse grid and restrict
      cout << "startamg" << endl;
      startamg(coarse,A,theta);
//      startamg(coarse,strhead,strnext,strprev,A,theta);
      cout << "makecoarsemx" << endl;
      makecoarsemx(coarse,coarse,A);
//      makecoarsemx(coarse,coarse,strhead,strnext,strprev,A);
      cout << "interpolate" << endl;
      interpolate(Acoarse,rcoarse,P,coarse,A,r);

// send to coarser level
      cout << "calling vcycle" << endl;
      xcoarse = new double[P.m];
      for (i = 0; i < P.m; i++)
         xcoarse[i] = 0.0;
//      outputAMGmx(delta2d,Acoarse);
      Vcycle(xcoarse,Acoarse,rcoarse,xcoarse,nu1,nu2,small);

// back from coarser level 
      cout << "back from vcycle" << endl;
      for (i = 0; i < P.n; i++)
         x[i] += getdotprod(P.row[i],P.nc[i],xcoarse);

// relax again
      cout << "gauss seidel" << endl;
      gaussseidel(x,A,b,x,nu2);

      delete [] rcoarse;
      delete [] xcoarse;
      clearsparse(Acoarse);
      clearsparse(P);
   }

   if (globnorm == 2)
      rnorm = getrelresidual(A,x,b);
   else
      rnorm = getrelresidualinf(A,x,b);

/*
   delta2d << A.n << " " << A.m << endl;
   delta2d << "x" << endl;
   for (i = 0; i < A.n; i++)
      delta2d << scientific << x[i] << " ";
   delta2d << endl;
*/

   delete [] coarse;
   delete [] r;
   return rnorm;
}
// used in amgsmall3
double Vcycle(double *x, SparseAMGList* &L, double *b, double *x0, int nu1, int nu2, 
              int small)
{
   int i, j, *coarse;
   SparseAMGList *currentlist;
   double rnorm, theta = 1.0; 
   double *rcoarse = NULL, *xcoarse = NULL, *r;
   int strhead, *strnext = NULL, *strprev = NULL;
   StrHeapStruct heap;
   
//   cout << ((*L).A).n << " " << ((*L).A).m << endl;

   for (i = 0; i < ((*L).A).n; i++)
      x[i] = x0[i];

   if (((*L).A).n <= small) 
   {
      if ((*L).B == NULL)
      {
         cout << ((*L).A).n << " " << ((*L).A).m << endl;
         ((*L).P).n = 0;
         ((*L).P).m = 0;
         (*L).B = matrix(((*L).A).n-1,((*L).A).m-1); 
         (*L).PLR = new int[((*L).A).n]; 
         (*L).PLC = new int[((*L).A).m];
      
         for (i = 0; i < ((*L).A).n; i++) 
            for (j = 0; j < ((*L).A).m; j++) 
               (*L).B[i][j] = 0.0;
         for (i = 0; i < ((*L).A).n; i++)
            for (j = 0; j < ((*L).A).nc[i]; j++) 
               (*L).B[i][(((*L).A).row[i][j]).roc] = (((*L).A).row[i][j]).val[0];
         gecp0((*L).B,(*L).PLR,(*L).PLC,(*L).B,((*L).A).n-1,((*L).A).m-1);
      }
      forwardbacksub0(x,b,(*L).B,(*L).PLR,(*L).PLC,((*L).A).n-1);
   }
   else 
   {
//      cout << "gauss seidel" << endl;
      gaussseidel(x,(*L).A,b,x,nu1);

// calculate residual
//      cout << "residual" << endl;
      r = new double[((*L).A).n];
      mxvecmult(r,(*L).A,x);
//      cout << "after mult" << endl;
      for (i = 0; i < ((*L).A).n; i++)
         r[i] = b[i]-r[i];

      if ((*L).next == NULL)
      {
         cout << ((*L).A).n << " " << ((*L).A).m << endl;
// make coarse grid and restrict
//         cout << "interpolate" << endl;
         coarse = new int[((*L).A).n];
         if (!globheap)
         {
            startamg(coarse,strhead,strnext,strprev,(*L).A,theta);
            makecoarsemx(coarse,coarse,strhead,strnext,strprev,(*L).A);
         }
         else
         {
            startamg(coarse,heap,(*L).A,theta);
            makecoarsemx(coarse,coarse,heap,(*L).A);
         }
         currentlist = new SparseAMGList;
         (*currentlist).next = NULL;
         (*currentlist).prev = L;
         (*currentlist).B = NULL;
         (*L).next = currentlist;
         interpolate((*currentlist).A, (*L).P, coarse, (*L).A); //(*currentlist).A is Ac,= = Pt A P
         delete [] coarse;

      }
      rcoarse = new double[((*L).P).m];
//      cout << "more residual" << endl;
      mxtvecmult(rcoarse,(*L).P,r); //Pt r
      delete [] r;

// send to coarser level
      xcoarse = new double[((*L).P).m];
      for (i = 0; i < ((*L).P).m; i++)
         xcoarse[i] = 0.0;
//      cout << "more vcycle" << endl;
      Vcycle(xcoarse,(*L).next,rcoarse,xcoarse,nu1,nu2,small);  // Solve for Ac ec = Pt r

// back from coarser level 
      for (i = 0; i < ((*L).P).n; i++)
         x[i] += getdotprod(((*L).P).row[i],((*L).P).nc[i],xcoarse); // correction x = x + P ec

// relax again
      gaussseidel(x,(*L).A,b,x,nu2);

      delete [] rcoarse;
      delete [] xcoarse;
   }

   if (globnorm == 2)
      rnorm = getrelresidual((*L).A,x,b);
   else
      rnorm = getrelresidualinf((*L).A,x,b);

/*
   delta2d << ((*L).A).n << " " << ((*L).A).m << endl;
   delta2d << "x" << endl;
   for (i = 0; i < ((*L).A).n; i++)
      delta2d << scientific << x[i] << " ";
   delta2d << endl;
*/

   return rnorm;
}

double getresidual(SparseAMGMx &A, double *x, double *b)
{
   int i;
   double temp, val = 0.0;

   for (i = 0; i < A.n; i++)
   {
      temp = getdotprod(A.row[i],A.nc[i],x);
//      if (fabs(b[i]-temp) > val)
//         val = fabs(b[i]-temp);
      val += (b[i]-temp)*(b[i]-temp);
   }
   val = sqrt(val);

   return val;
}

double getresidualinf(SparseAMGMx &A, double *x, double *b)
{
   int i;
   double temp, val = 0.0;

   for (i = 0; i < A.n; i++)
   {
      temp = getdotprod(A.row[i],A.nc[i],x);
      if (fabs(b[i]-temp) > val)
         val = fabs(b[i]-temp);
   }

   return val;
}

double getrelresidual(SparseAMGMx &A, double *x, double *b)
{
   int i;
   double themax = 0.0;

   for (i = 0; i < A.m; i++)
//      if (fabs(b[i]) > themax)
//         themax = fabs(b[i]);
      themax += b[i]*b[i];
   themax = sqrt(themax);

   if (themax != 0.0)
      return getresidual(A,x,b)/themax;
   else
      return getresidual(A,x,b);
}

double getrelresidualinf(SparseAMGMx &A, double *x, double *b)
{
   int i;
   double themax = 0.0;

   for (i = 0; i < A.m; i++)
      if (fabs(b[i]) > themax)
         themax = fabs(b[i]);

   if (themax != 0.0)
      return getresidualinf(A,x,b)/themax;
   else
      return getresidualinf(A,x,b);
}

void Vcycleloop(double *x, SparseAMGMx &A, double *b, double *x0, int nu1, int nu2, 
                int small, int numsteps, double tol)
{
   int i, step = 0;
   double residual;

   cout << "in loop" << endl;
   for (i = 0; i < A.n; i++)
      x[i] = x0[i];
   cout << "getting residual" << endl;
   if (globnorm == 2)
      residual = getrelresidual(A,x,b);
   else
      residual = getrelresidualinf(A,x,b);
//   cout.precision(16);
//   cout << "step = " << step << " and residual = " << scientific << residual << endl;
   cout << "step = " << step << " and residual = " << residual << endl;
   for (step = 1; step <= numsteps && residual > tol; step++)
   {
      residual = Vcycle(x,A,b,x,nu1,nu2,small);
//      cout << "step = " << step << " and residual = " << scientific << residual << endl;
      cout << "step = " << step << " and residual = " << residual << endl;
   }

   cout << "AMG took " << step-1 << " steps to converge." << endl;
}

void Vcycleloop(double *x, SparseAMGList* &L, double *b, double *x0, int nu1, int nu2, 
                int small, int numsteps, double tol)
{
   int i, step = 0;
   double residual;

   for (i = 0; i < ((*L).A).n; i++)
      x[i] = x0[i];
   if (globnorm == 2)
      residual = getrelresidual((*L).A,x,b);
   else
      residual = getrelresidualinf((*L).A,x,b);
//   cout.precision(16);
//   cout << "step = " << step << " and residual = " << scientific << residual << endl;
   cout << "step = " << step << " and residual = " << residual << endl;
   for (step = 1; step <= numsteps && residual > tol; step++)
   {
      residual = Vcycle(x,L,b,x,nu1,nu2,small);
      if (step%10 == 0)
         cout << "step = " << step << " and residual = " << residual << endl;
   }

//   cout << "AMG took " << step-1 << " steps to converge with residual " << scientific 
//        << residual << endl;
   cout << "AMG took " << step-1 << " steps to converge with residual " 
        << residual << endl;
}

void AMGsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, int nu1, 
              int nu2, int small, int numsteps, double tol, double ***S, PBData &pb)
{
// algebraic multigrid that does not remove beforehand identity elements in the matrix
   int i, j, r, m, n, count, N = 1, tindex[grid.dim], rindex[grid.dim];
   for (i = 0; i < grid.dim; i++)
      N *= grid.nx[i]+1;
   double ehere, ethere, value, amgx[N], amgb[N];
   SparseElt *temp = NULL, *current; 
   SparseElt2 *current2;

   cout << "here 1" << endl;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      amgx[i] = evalarray(x,tindex);
      amgb[i] = evalarray(b,tindex);
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "here 2" << endl;
   SparseAMGMx amgA;
   amgA.n = N;
   amgA.m = N;
   amgA.row = new SparseAMGElt*[amgA.n];
   amgA.nc = new int[amgA.n];
   amgA.col = new SparseAMGElt*[amgA.m];
   amgA.nr = new int[amgA.m];
   SparseElt **amgAt = new SparseElt*[amgA.m];
   for (i = 0; i < amgA.m; i++)
      amgAt[i] = NULL;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      temp = NULL;
      if (evalarray(A,tindex) == NULL)
      {
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
         value = 0.0;
         for (m = 0; m < grid.dim; m++)
            value += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
         sparseorder(i,temp,value);
         count = 1;

         for (m = 0; m < grid.dim; m++)
            for (n = -1; n <= 1; n += 2)
            {
               for (r = 0; r < grid.dim; r++) 
                  rindex[r] = tindex[r];
               rindex[m] = tindex[m]+n;
               j = sub2ind(rindex,grid.nx,grid.dim);
               if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
               {
                  sparseorder(j,temp,-ehere/(grid.dx[m]*grid.dx[m]));
                  count++;
               }
            }
      }
      else
      {
         count = 0;
         for (current2 = evalarray(A,tindex); current2 != NULL; 
              current2 = (*current2).next)
         {
            j = sub2ind((*current2).cindex,grid.nx,grid.dim);
            sparseorder(j,temp,(*current2).val);
            count++;
         }
      }

      amgA.nc[i] = count;
      amgA.row[i] = new SparseAMGElt[amgA.nc[i]];
      count = 0;
      for (current = temp; current != NULL; 
           temp = current, current = (*current).next, delete temp)
      {
         (amgA.row[i][count]).roc = (*current).cindex;
         (amgA.row[i][count]).val = new double[1];
         (amgA.row[i][count]).val[0] = (*current).val;
         (amgA.row[i][count]).strong = new char[1];
         (amgA.row[i][count]).strong[0] = 0;
         sparseorder((amgA.row[i][count]).roc,i,amgAt,(amgA.row[i][count]).val[0]);
         count++;
      }
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   for (i = 0; i < amgA.m; i++)
   {
      count = 0;
      for (current = amgAt[i]; current != NULL; current = (*current).next)
         count++;
      amgA.nr[i] = count;
      amgA.col[i] = new SparseAMGElt[amgA.nr[i]];
      count = 0;
      for (current = amgAt[i]; current != NULL; 
           temp = current, current = (*current).next, delete temp)
      {
         (amgA.col[i][count]).roc = (*current).cindex;
         j = (amgA.col[i][count]).roc;
         for (r = 0; r < amgA.nc[j] && (amgA.row[j][r]).roc != i; r++);
         (amgA.col[i][count]).val = (amgA.row[j][r]).val;
         (amgA.col[i][count]).strong = (amgA.row[j][r]).strong;
         count++;
      }
   }
   delete[] amgAt;

   cout << "here 3" << endl;
   Vcycleloop(amgx,amgA,amgb,amgx,nu1,nu2,small,numsteps,tol);
   clearsparse(amgA);

   cout << "here 4" << endl;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      setvalarray(x,tindex,amgx[i]);
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cout << "here 5" << endl;
}

void AMGsmall2(double ***x, SparseElt2**** &A, double ***b, GridData &grid, int nu1, 
               int nu2, int small, int numsteps, double tol, double ***S, PBData &pb)
{
// algebraic multigrid that removes beforehand identity elements in the matrix
// recomputes matrices everytime
   int i, j, r, m, n, count, N = 1, tindex[grid.dim], rindex[grid.dim];
   for (i = 0; i < grid.dim; i++)
      N *= grid.nx[i]+1;
   double ehere, ethere, value, *amgx = new double[N], *amgb = new double[N];
   int *position = new int[N];
   char *nontriv = new char[N];
   SparseElt *temp = NULL, *current; 
   SparseElt2 *current2;

   cout << "here 1" << endl;
   N = 0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      count = 0;
      for (current2 = evalarray(A,tindex); current2 != NULL; 
           current2 = (*current2).next)
         if (sub2ind((*current2).cindex,grid.nx,grid.dim) != i)
            count++;
         else
            value = (*current2).val;
      if (evalarray(A,tindex) != NULL && count == 0)
      {
         setvalarray(x,tindex,evalarray(b,tindex)/value);
         nontriv[i] = 0;
         position[i] = -1;
      }
      else
      {
         amgx[N] = evalarray(x,tindex);
         amgb[N] = evalarray(b,tindex);
         nontriv[i] = 1;
         position[i] = N;
         N++;
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "N = " << N << " " << (grid.nx[0]+1)*(grid.nx[1]+1)*(grid.nx[2]+1)-N << endl;
   SparseAMGMx amgA;
   amgA.n = N;
   amgA.m = N;
   amgA.row = new SparseAMGElt*[amgA.n];
   amgA.nc = new int[amgA.n];
   amgA.col = new SparseAMGElt*[amgA.m];
   amgA.nr = new int[amgA.m];
   SparseElt **amgAt = new SparseElt*[amgA.m];
   for (i = 0; i < amgA.m; i++)
      amgAt[i] = NULL;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      if (nontriv[i])
      {
         temp = NULL;
         if (evalarray(A,tindex) == NULL)
         {
            if (evalarray(S,tindex) < 0.0)
               ehere = pb.epsilonm;
            else
               ehere = pb.epsilonp;
            value = 0.0;
            for (m = 0; m < grid.dim; m++)
               value += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
            sparseorder(position[i],temp,value);
            count = 1;

            for (m = 0; m < grid.dim; m++)
               for (n = -1; n <= 1; n += 2)
               {
                  for (r = 0; r < grid.dim; r++) 
                     rindex[r] = tindex[r];
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                  {
                     j = sub2ind(rindex,grid.nx,grid.dim);
                     value = -ehere/(grid.dx[m]*grid.dx[m]);
                     if (nontriv[j])
                     {
                        sparseorder(position[j],temp,value);
                        count++;
                     }
                     else
                        amgb[position[i]] -= value*evalarray(x,rindex);
                  }
               }
         }
         else
         {
            count = 0;
            for (current2 = evalarray(A,tindex); current2 != NULL; 
                 current2 = (*current2).next)
            {
               j = sub2ind((*current2).cindex,grid.nx,grid.dim);
               if (nontriv[j])
               {
                  sparseorder(position[j],temp,(*current2).val);
                  count++;
               }
               else
                  amgb[position[i]] -= (*current2).val*evalarray(x,(*current2).cindex);
            }
         }

         amgA.nc[position[i]] = count;
         amgA.row[position[i]] = new SparseAMGElt[amgA.nc[position[i]]];
         count = 0;
         for (current = temp; current != NULL; 
              temp = current, current = (*current).next, delete temp)
         {
            (amgA.row[position[i]][count]).roc = (*current).cindex;
            (amgA.row[position[i]][count]).val = new double[1];
            (amgA.row[position[i]][count]).val[0] = (*current).val;
            (amgA.row[position[i]][count]).strong = new char[1];
            (amgA.row[position[i]][count]).strong[0] = 0;
            sparseorder((amgA.row[position[i]][count]).roc,position[i],amgAt,
                        (amgA.row[position[i]][count]).val[0]);
            count++;
         }
      }
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "here 3" << endl;
   for (i = 0; i < amgA.m; i++)
   {
      count = 0;
      for (current = amgAt[i]; current != NULL; current = (*current).next)
         count++;
      amgA.nr[i] = count;
      amgA.col[i] = new SparseAMGElt[amgA.nr[i]];
      count = 0;
      for (current = amgAt[i]; current != NULL; 
           temp = current, current = (*current).next, delete temp)
      {
         (amgA.col[i][count]).roc = (*current).cindex;
         j = (amgA.col[i][count]).roc;
         for (r = 0; r < amgA.nc[j] && (amgA.row[j][r]).roc != i; r++);
         (amgA.col[i][count]).val = (amgA.row[j][r]).val;
         (amgA.col[i][count]).strong = (amgA.row[j][r]).strong;
         count++;
      }
   }
   delete[] amgAt;

   cout << "here 3" << endl;
//   outputAMGmx(delta2d,amgA);
   Vcycleloop(amgx,amgA,amgb,amgx,nu1,nu2,small,numsteps,tol);
   clearsparse(amgA);

   cout << "here 4" << endl;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      if (nontriv[i])
         setvalarray(x,tindex,amgx[position[i]]);
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cout << "here 5" << endl;

   delete [] amgx; 
   delete [] amgb;
   delete [] position;
   delete [] nontriv;
}
// used ing getRHS
// nu1 = 2, nu2 = 0, small = grid.nx[0], numsteps = grid.nx[0]^3
// , initilized as 0
void AMGsmall3(double ***x, SparseElt2**** &A, double ***b, GridData &grid, int nu1, 
               int nu2, int small, int numsteps, double tol, double ***S, PBData &pb)
{
// algebraic multigrid that removes beforehand identity elements in the matrix
// stores matrices for use in subsequent interations
   int i, j, r, m, n, count, N = 1, tindex[grid.dim], rindex[grid.dim];
   for (i = 0; i < grid.dim; i++)
      N *= grid.nx[i]+1;
   double ehere, ethere, value, *amgx = new double[N], *amgb = new double[N];
   int *position = new int[N];
   char *nontriv = new char[N];
   SparseElt *temp = NULL, *current; 
   SparseElt2 *current2;

   N = 0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      count = 0;
      for (current2 = evalarray(A,tindex); current2 != NULL; 
           current2 = (*current2).next)
         if (sub2ind((*current2).cindex,grid.nx,grid.dim) != i)
            count++; // count number of nonzeros in a row, excluding diaganal
         else
            value = (*current2).val; // get diagonal value
      if (evalarray(A,tindex) != NULL && count == 0) // boundary points
      {
         setvalarray(x,tindex,evalarray(b,tindex)/value); //initial guess is b
         nontriv[i] = 0;
         position[i] = -1;
      }
      else // for interface point and interior points
      {
         amgx[N] = evalarray(x,tindex);
         amgb[N] = evalarray(b,tindex);
         nontriv[i] = 1;
         position[i] = N;
         N++;
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   SparseAMGList *L, *currentlist;
   L = new SparseAMGList;
   ((*L).next) = NULL;
   ((*L).prev) = NULL;
   ((*L).A).n = N; // N is the number of nontrivial element
   ((*L).A).m = N;
   ((*L).A).row = new SparseAMGElt*[((*L).A).n];
   ((*L).A).nc = new int[((*L).A).n];
   ((*L).A).col = new SparseAMGElt*[((*L).A).m];
   ((*L).A).nr = new int[((*L).A).m];
   SparseElt **amgAt = new SparseElt*[((*L).A).m]; // N by N matrix of sparse element
   for (i = 0; i < ((*L).A).m; i++)
      amgAt[i] = NULL;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      if (nontriv[i]) // interface or interior
      {
         temp = NULL;
         if (evalarray(A,tindex) == NULL)// interior points
         {
            if (evalarray(S,tindex) < 0.0)
               ehere = pb.epsilonm;
            else
               ehere = pb.epsilonp;
            value = 0.0;
            for (m = 0; m < grid.dim; m++)
               value += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
            sparseorder(position[i],temp,value);
            count = 1;

            for (m = 0; m < grid.dim; m++)
               for (n = -1; n <= 1; n += 2)
               {
                  for (r = 0; r < grid.dim; r++) 
                     rindex[r] = tindex[r];
                  rindex[m] = tindex[m]+n;
                  j = sub2ind(rindex,grid.nx,grid.dim);
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                  {
                     value = -ehere/(grid.dx[m]*grid.dx[m]);
                     if (nontriv[j]) // j is not boundary
                     {
                        sparseorder(position[j],temp,value);
                        count++;
                     }
                     else
                        amgb[position[i]] -= value*evalarray(x,rindex); //x is zero here
                  }
               }
         }
         else
         {
            count = 0;
            for (current2 = evalarray(A,tindex); current2 != NULL; 
                 current2 = (*current2).next)
            {
               j = sub2ind((*current2).cindex,grid.nx,grid.dim);
               if (nontriv[j]) //nontriv[j] == 1, j is interior
               {
                  sparseorder(position[j],temp,(*current2).val);
                  count++;
               }
               else
                  amgb[position[i]] -= (*current2).val*evalarray(x,(*current2).cindex);
            }
         }

         ((*L).A).nc[position[i]] = count;
         ((*L).A).row[position[i]] = new SparseAMGElt[((*L).A).nc[position[i]]];
         count = 0;
         for (current = temp; current != NULL; 
              temp = current, current = (*current).next, delete temp)
         {
            (((*L).A).row[position[i]][count]).roc = (*current).cindex;
            (((*L).A).row[position[i]][count]).val = new double[1];
            (((*L).A).row[position[i]][count]).val[0] = (*current).val;
            (((*L).A).row[position[i]][count]).strong = new char[1];
            (((*L).A).row[position[i]][count]).strong[0] = 0;
            sparseorder((((*L).A).row[position[i]][count]).roc, position[i],amgAt,
                        (((*L).A).row[position[i]][count]).val[0]); //
            count++;
         }
      }
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   for (i = 0; i < ((*L).A).m; i++)
   {
      count = 0;
      for (current = amgAt[i]; current != NULL; current = (*current).next)
         count++;
      ((*L).A).nr[i] = count;
      ((*L).A).col[i] = new SparseAMGElt[((*L).A).nr[i]];
      count = 0;
      for (current = amgAt[i]; current != NULL; 
           temp = current, current = (*current).next, delete temp)
      {
         (((*L).A).col[i][count]).roc = (*current).cindex;
         j = (((*L).A).col[i][count]).roc;
         for (r = 0; r < ((*L).A).nc[j] && (((*L).A).row[j][r]).roc != i; r++);
         (((*L).A).col[i][count]).val = (((*L).A).row[j][r]).val;
         (((*L).A).col[i][count]).strong = (((*L).A).row[j][r]).strong;
         count++;
      }
   }
   
   delete[] amgAt;
   // OutputAmgMxByRow("amgrow.txt",(*L).A);
   // OutputVector("amgb.txt",amgb, N);
   // OutputAmgMxByCol("amgcol.txt",(*L).A);
//   outputAMGmx(delta2d,(*L).A);
   Vcycleloop(amgx,L,amgb,amgx,nu1,nu2,small,numsteps,tol);
   for ( ; L != NULL; 
        currentlist = L, L = (*L).next, delete currentlist) 
   {
      clearsparse((*L).A);
      clearsparse((*L).P);
      if ((*L).next == NULL)
      {
         free_matrix((*L).B,((*L).A).n-1,((*L).A).n-1);
         delete [] (*L).PLR;
         delete [] (*L).PLC;
      }
   }

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      if (nontriv[i])
         setvalarray(x,tindex,amgx[position[i]]);
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   delete [] amgx; 
   delete [] amgb;
   delete [] position;
   delete [] nontriv;
}
// with a
void AMGsmall3(double ***x, SparseElt2**** &A, double ***b, GridData &grid, int nu1, 
               int nu2, int small, int numsteps, double tol, double ***a, double ***S, 
               PBData &pb)
{
// algebraic multigrid that removes beforehand identity elements in the matrix
// stores matrices for use in subsequent interations
   int i, j, r, m, n, count, N = 1, tindex[grid.dim], rindex[grid.dim];
   for (i = 0; i < grid.dim; i++)
      N *= grid.nx[i]+1;
   double ehere, ethere, value, *amgx = new double[N], *amgb = new double[N];
   int *position = new int[N];
   char *nontriv = new char[N];
   SparseElt *temp = NULL, *current; 
   SparseElt2 *current2;

   N = 0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      count = 0;
      for (current2 = evalarray(A,tindex); current2 != NULL; 
           current2 = (*current2).next)
         if (sub2ind((*current2).cindex,grid.nx,grid.dim) != i)
            count++;
         else
            value = (*current2).val;
      if (evalarray(A,tindex) != NULL && count == 0)
      {
         setvalarray(x,tindex,evalarray(b,tindex)/value);
         nontriv[i] = 0;
         position[i] = -1;
      }
      else
      {
         amgx[N] = evalarray(x,tindex);
         amgb[N] = evalarray(b,tindex);
         nontriv[i] = 1;
         position[i] = N;
         N++;
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   SparseAMGList *L, *currentlist;
   L = new SparseAMGList;
   ((*L).next) = NULL;
   ((*L).prev) = NULL;
   ((*L).A).n = N;
   ((*L).A).m = N;
   ((*L).A).row = new SparseAMGElt*[((*L).A).n];
   ((*L).A).nc = new int[((*L).A).n];
   ((*L).A).col = new SparseAMGElt*[((*L).A).m];
   ((*L).A).nr = new int[((*L).A).m];
   SparseElt **amgAt = new SparseElt*[((*L).A).m];
   for (i = 0; i < ((*L).A).m; i++)
      amgAt[i] = NULL;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      if (nontriv[i])
      {
         temp = NULL;
         if (evalarray(A,tindex) == NULL)
         {
            if (evalarray(S,tindex) < 0.0)
               ehere = pb.epsilonm;
            else
               ehere = pb.epsilonp;
            value = evalarray(a,tindex);
            for (m = 0; m < grid.dim; m++)
               value += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
            sparseorder(position[i],temp,value);
            count = 1;

            for (m = 0; m < grid.dim; m++)
               for (n = -1; n <= 1; n += 2)
               {
                  for (r = 0; r < grid.dim; r++) 
                     rindex[r] = tindex[r];
                  rindex[m] = tindex[m]+n;
                  j = sub2ind(rindex,grid.nx,grid.dim);
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                  {
                     value = -ehere/(grid.dx[m]*grid.dx[m]);
                     if (nontriv[j])
                     {
                        sparseorder(position[j],temp,value);
                        count++;
                     }
                     else
                        amgb[position[i]] -= value*evalarray(x,rindex);
                  }
               }
         }
         else
         {
            count = 0;
            for (current2 = evalarray(A,tindex); current2 != NULL; 
                 current2 = (*current2).next)
            {
               j = sub2ind((*current2).cindex,grid.nx,grid.dim);
               if (nontriv[j])
               {
                  sparseorder(position[j],temp,(*current2).val);
                  count++;
               }
               else
                  amgb[position[i]] -= (*current2).val*evalarray(x,(*current2).cindex);
            }
         }

         ((*L).A).nc[position[i]] = count;
         ((*L).A).row[position[i]] = new SparseAMGElt[((*L).A).nc[position[i]]];
         count = 0;
         for (current = temp; current != NULL; 
              temp = current, current = (*current).next, delete temp)
         {
            (((*L).A).row[position[i]][count]).roc = (*current).cindex;
            (((*L).A).row[position[i]][count]).val = new double[1];
            (((*L).A).row[position[i]][count]).val[0] = (*current).val;
            (((*L).A).row[position[i]][count]).strong = new char[1];
            (((*L).A).row[position[i]][count]).strong[0] = 0;
            sparseorder((((*L).A).row[position[i]][count]).roc,position[i],amgAt,
                        (((*L).A).row[position[i]][count]).val[0]);
            count++;
         }
      }
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   for (i = 0; i < ((*L).A).m; i++)
   {
      count = 0;
      for (current = amgAt[i]; current != NULL; current = (*current).next)
         count++;
      ((*L).A).nr[i] = count;
      ((*L).A).col[i] = new SparseAMGElt[((*L).A).nr[i]];
      count = 0;
      for (current = amgAt[i]; current != NULL; 
           temp = current, current = (*current).next, delete temp)
      {
         (((*L).A).col[i][count]).roc = (*current).cindex;
         j = (((*L).A).col[i][count]).roc;
         for (r = 0; r < ((*L).A).nc[j] && (((*L).A).row[j][r]).roc != i; r++);
         (((*L).A).col[i][count]).val = (((*L).A).row[j][r]).val;
         (((*L).A).col[i][count]).strong = (((*L).A).row[j][r]).strong;
         count++;
      }
   }
   delete[] amgAt;

//   outputAMGmx(delta2d,(*L).A);
   Vcycleloop(amgx,L,amgb,amgx,nu1,nu2,small,numsteps,tol);
   for ( ; L != NULL; 
        currentlist = L, L = (*L).next, delete currentlist) 
   {
      clearsparse((*L).A);
      clearsparse((*L).P);
      if ((*L).next == NULL)
      {
         free_matrix((*L).B,((*L).A).n-1,((*L).A).n-1);
         delete [] (*L).PLR;
         delete [] (*L).PLC;
      }
   }

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      i = sub2ind(tindex,grid.nx,grid.dim);
      if (nontriv[i])
         setvalarray(x,tindex,amgx[position[i]]);
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   delete [] amgx; 
   delete [] amgb;
   delete [] position;
   delete [] nontriv;
}

void checkmx(SparseAMGMx A)
{
   int i, j, r, s;

   for (i = 0; i < A.n; i++)
      for (r = 0; r < A.nc[i]-1; r++)
         if ((A.row[i][r]).roc >= (A.row[i][r+1]).roc)
         {
            cout << "row order is wrong" << endl;
            exit(1);
         }

   for (i = 0; i < A.m; i++)
      for (r = 0; r < A.nr[i]-1; r++)
         if ((A.col[i][r]).roc >= (A.col[i][r+1]).roc)
         {
            cout << "column order is wrong: " << i << " " << r << " "
                 << (A.col[i][r]).roc << " " << (A.col[i][r+1]).roc << endl;
            exit(1);
         }

   for (i = 0; i < A.n; i++)
      for (r = 0; r < A.nc[i]; r++)
      {
         j = (A.row[i][r]).roc;
         for (s = 0; s < A.nr[j] && (A.col[j][s]).roc != i; s++);
         if (s >= A.nr[j])
         {
            cout << "row element does not exist in column: " << i << " " << r << " " 
                 << (A.row[i][r]).roc << " " << j << endl;
            exit(1);
         }
      }

   for (i = 0; i < A.m; i++)
      for (r = 0; r < A.nr[i]; r++)
      {
         j = (A.col[i][r]).roc;
         for (s = 0; s < A.nc[j] && (A.row[j][s]).roc != i; s++);
         if (s >= A.nc[j])
         {
            cout << "column element does not exist in row" << endl;
            exit(1);
         }
      }

   for (i = 0; i < A.n; i++)
      for (r = 0; r < A.nc[i]; r++)
      {
         j = (A.row[i][r]).roc;
         for (s = 0; s < A.nr[j] && (A.col[j][s]).roc != i; s++);
         if ((A.row[i][r]).val != (A.col[j][s]).val)
         {
            cout << "row and column values do not correspond" << endl;
            exit(1);
         }
      }

   for (i = 0; i < A.n; i++)
      for (r = 0; r < A.nc[i]; r++)
      {
         j = (A.row[i][r]).roc;
         for (s = 0; s < A.nr[j] && (A.col[j][s]).roc != i; s++);
         if ((A.row[i][r]).strong != (A.col[j][s]).strong)
         {
            cout << "row and column strengths do not correspond" << endl;
            exit(1);
         }
      }

   cout << "Matrix checks out" << endl;
}

void peekmx(SparseAMGMx A)
{
   int i, r, count;
   double value, tol = 1.0e-14;

   count = 0;
   for (i = 0; i < A.n; i++)
      for (r = 0; r < A.nc[i]; r++) 
         if ((A.row[i][r]).roc != i && (A.row[i][r]).val[0] > 0.0)
            count++;

   cout << "Number of off-diagonal elements of wrong sign = " << count << endl;

   count = 0;
   for (i = 0; i < A.n; i++)
      for (r = 0; r < A.nc[i]; r++) 
         if ((A.row[i][r]).roc == i && (A.row[i][r]).val[0] <= 0.0)
            count++;

   cout << "Number of diagonal elements of wrong sign = " << count << endl;

   count = 0;
   for (i = 0; i < A.n; i++)
   {
      value = 0.0;
      for (r = 0; r < A.nc[i]; r++) 
         if ((A.row[i][r]).roc != i)
            value -= fabs((A.row[i][r]).val[0]);
         else if ((A.row[i][r]).roc == i)
            value += (A.row[i][r]).val[0];
      if (value < -tol)
         count++;
   }

   cout << "Number of rows not diagonally dominant = " << count << endl;

   count = 0;
   for (i = 0; i < A.n; i++)
   {
      value = 0.0;
      for (r = 0; r < A.nc[i]; r++) 
         if ((A.row[i][r]).roc != i)
            value += fabs((A.row[i][r]).val[0]);
      if (value < tol)
         count++;
   }

   cout << "Number of identity rows = " << count << endl;
}

void getpoisson(SparseAMGMx &A, double *b, int *nx, double *dx, int thedim)
{
   int i, r, s, t, row, count, temp; 
   int tindex[thedim], rindex[thedim], col[2*thedim+1];
   
   cout << "getpoisson" << endl;
   A.n = 1.0;
   for (r = 0; r < thedim; r++)
      A.n *= nx[r]+1;
   A.m = A.n;
   A.row = new SparseAMGElt*[A.n];
   A.nc = new int[A.n];
   A.col = new SparseAMGElt*[A.m];
   A.nr = new int[A.m];
   for (i = 0; i < A.n; i++)
   {
      A.nc[i] = 0;
      A.row[i] = new SparseAMGElt[2*thedim+1];
   }
   for (i = 0; i < A.m; i++)
   {
      A.nr[i] = 0;
      A.col[i] = new SparseAMGElt[2*thedim+1];
   }
   
   for (i = 0; i < thedim; i++)
      tindex[i] = 0;
   while (tindex[0] <= nx[0])
   {
      row = sub2ind(tindex,nx,thedim);
      b[row] = 1.0;

      for (r = 0; r < thedim; r++)
         rindex[r] = tindex[r];
      col[0] = row;
      count = 1;
      for (r = 0; r < thedim; r++)
      {
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = tindex[r]+s;
            if (rindex[r] >= 0 && rindex[r] <= nx[r])
            {
               col[count] = sub2ind(rindex,nx,thedim);
               count++;
               for (t = count-1; t >= 1; t--)
                  if (col[t] < col[t-1])
                  {
                     temp = col[t];
                     col[t] = col[t-1];
                     col[t-1] = temp;
                  }
            }
         }
         rindex[r] = tindex[r];
      }

      for (t = 0; t < count; t++)
      {
         (A.row[row][t]).roc = col[t];
         (A.row[row][t]).strong = new char[1];
         (A.row[row][t]).strong[0] = 0;
         (A.row[row][t]).val = new double[1];
         if (col[t] != row)
            (A.row[row][t]).val[0] = -1.0/(dx[0]*dx[0]);
         else
            (A.row[row][t]).val[0] = 2.0*thedim/(dx[0]*dx[0]);

         (A.col[col[t]][A.nr[col[t]]]).roc = row;
         (A.col[col[t]][A.nr[col[t]]]).strong = (A.row[row][t]).strong;
         (A.col[col[t]][A.nr[col[t]]]).val = (A.row[row][t]).val;

         (A.nc[row])++;
         (A.nr[col[t]])++;
      }

      (tindex[thedim-1])++;
      for (i = thedim-1; i > 0 && tindex[i] > nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}

char closertohead(StrHeapStruct &heap, int r, int s)
{
   if (heap.value[heap.index[r]] > heap.value[heap.index[s]])
      return 1;
   else if (heap.value[heap.index[r]] < heap.value[heap.index[s]])
      return 0;
   else if (heap.index[r] < heap.index[s])
      return 1;
   else
      return 0;
}

void fixheap(StrHeapStruct &heap, int num)
{
   fixheapdown(heap,num);
   fixheapup(heap,num);
}

void fixheapdown(StrHeapStruct &heap, int num)
{
   int r;

   for (r = num; 2*r+1 < heap.num && 2*r+2 < heap.num && 
                 (closertohead(heap,2*r+1,r) || closertohead(heap,2*r+2,r)); )
      if (closertohead(heap,2*r+1,2*r+2))
      {
         switchheapelt(heap,r,2*r+1);
         r = 2*r+1;
      }
      else
      {
         switchheapelt(heap,r,2*r+2);
         r = 2*r+2;
      }
   if (2*r+1 < heap.num && 2*r+2 >= heap.num && closertohead(heap,2*r+1,r))
      switchheapelt(heap,r,2*r+1);
}

void fixheapup(StrHeapStruct &heap, int num)
{
   int r;

   for (r = num; r > 0 && closertohead(heap,r,(r-1)/2); r = (r-1)/2)
      switchheapelt(heap,r,(r-1)/2);
}

void addtoheap(StrHeapStruct &heap, int index)
{
   heap.index[heap.num] = index;
   heap.loc[heap.index[heap.num]] = heap.num;
   (heap.num)++;
   fixheapup(heap,heap.num-1);
}

void popfromheap(StrHeapStruct &heap, int num)
{
   int r;

   for (r = num; 2*r+1 < heap.num && 2*r+2 < heap.num; )
      if (closertohead(heap,2*r+1,2*r+2))
      {
         switchheapelt(heap,r,2*r+1);
         r = 2*r+1;
      }
      else
      {
         switchheapelt(heap,r,2*r+2);
         r = 2*r+2;
      }
   if (2*r+1 < heap.num)
   {
      switchheapelt(heap,r,2*r+1);
      r = 2*r+1;
   }
  
   switchheapelt(heap,r,heap.num-1);
   (heap.num)--;
   if (r < heap.num)
      fixheapup(heap,r);
}

void switchheapelt(StrHeapStruct &heap, int r, int s)
{
   int temp;

   temp = heap.index[r];
   heap.index[r] = heap.index[s];
   heap.index[s] = temp;

   heap.loc[heap.index[r]] = r;
   heap.loc[heap.index[s]] = s;
}

void sortfromheap(StrHeapStruct &heap)
{
   while (heap.num > 0)
      popfromheap(heap,0);
}

char checkheap(StrHeapStruct &heap)
{
   int r;

   for (r = 0; r < heap.num; r++)
   {
      if (2*r+1 < heap.num && !closertohead(heap,r,2*r+1))
      {
         cout << "Error at " << r << " with child " << 2*r+1 << endl;
         return 0;
      }
      if (2*r+2 < heap.num && !closertohead(heap,r,2*r+2))
      {
         cout << "Error at " << r << " with child " << 2*r+2 << endl;
         return 0;
      }
   }
   return 1;
}


void outputAMGmx(ofstream& outfile, SparseAMGMx A)
{
   int i, r;

   outfile.precision(16);
   for (i = 0; i < A.n; i++)
      for (r = 0; r < A.nc[i]; r++)
         outfile << i << " " << (A.row[i][r]).roc << " " << scientific 
                 << (A.row[i][r]).val[0] << endl;
}


void OutputAmgMxByCol(string filename, SparseAMGMx A)
{
   int i, r;
   ofstream outfile(filename);
   outfile.precision(16);
   for (i = 0; i < A.m; i++)
      for (r = 0; r < A.nr[i]; r++)
         outfile << i << " " << (A.col[i][r]).roc << " " << scientific 
                 << (A.col[i][r]).val[0] << endl;
}

void OutputAmgMxByRow(string filename, SparseAMGMx A)
{
   int i, r;
   ofstream outfile(filename);
   outfile.precision(16);
   for (i = 0; i < A.n; i++)
      for (r = 0; r < A.nc[i]; r++)
         outfile << i << " " << (A.row[i][r]).roc << " " << scientific 
                 << (A.row[i][r]).val[0] << endl;
}