#include "sparse.h"
#include <iostream>
using namespace std;

SparseElt2 ****sparseelt2ptrmatrix(int row, int col, int fr)
{
   int i, j, k;
   SparseElt2 ****m;
   int err = 0;

   m = new SparseElt2***[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new SparseElt2**[col+1];
      if (!m[i] && !err)
         err = 1;
      for (j = 0; j <= col; j++)
      {
         m[i][j] = new SparseElt2*[fr+1];
         if (!m[i][j] && !err)
            err = 1;
         for (k = 0; k <= fr; k++)
            m[i][j][k] = NULL;
      }
   }
   // if (err)
   //    cout << "matrix stuff didn`t work" << endl;
   return m;
}



//rindex is row index
//cindex is column index
//always append
void sparse2(int *rindex, int *cindex, SparseElt2**** &A, double value, GridData &grid)
{
   int r;
   SparseElt2 *current;
   //current = A[rindex], if current is not NULL and current.cindex is not the same as cindex, go to next
   for (current = evalarray(A,rindex); current != NULL && 
        sub2ind((*current).cindex,grid.nx,grid.dim) != sub2ind(cindex,grid.nx,grid.dim); 
        current = (*current).next);
   if (current == NULL)//if current is NULL, A[rindex].cindex = cindex, A[rindex].next = A[rindex]
   {
      SparseElt2 *temp;
      temp = new SparseElt2;
      (*temp).cindex = new int[grid.dim];
      for (r = 0; r < grid.dim; r++) 
         (*temp).cindex[r] = cindex[r];
      (*temp).val = value;
      (*temp).next = evalarray(A,rindex);
      setvalarray(A,rindex,temp);
   }
   else
      (*current).val += value;
}

/*
A is a linked list of sparse element
if cindex is in A, set A[cindex] += value
otherwise A[cindex]= value, append at the front
*/ 
void sparse2(int *cindex, SparseElt2* &A, double value, GridData &grid)
{
   int r;
   SparseElt2 *current;

   for (current = A; current != NULL && 
        sub2ind((*current).cindex,grid.nx,grid.dim) != sub2ind(cindex,grid.nx,grid.dim); 
        current = (*current).next);
   if (current == NULL)
   {
      SparseElt2 *temp;
      temp = new SparseElt2;
      (*temp).cindex = new int[grid.dim];
      for (r = 0; r < grid.dim; r++) 
         (*temp).cindex[r] = cindex[r];
      (*temp).val = value;
      (*temp).next = A;
      A = temp;
   }
   else
      (*current).val += value;
}

void sparse2(int *index, SparseElt2* &R, SparseElt2**** &A, GridData &grid)
{
   SparseElt2 *current;

   for (current = R; current != NULL; current = (*current).next)
      sparse2(index,(*current).cindex,A,(*current).val,grid);
}

void sparseorder(int *rindex, int *cindex, SparseElt2**** &A, double value, 
                 GridData &grid)
{
   int r;
   SparseElt2 *current, *prev;

   for (prev = NULL,current = evalarray(A,rindex); 
        current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) <=
                           sub2ind(cindex,grid.nx,grid.dim); 
        prev = current,current = (*current).next);
   if (prev == NULL || sub2ind((*prev).cindex,grid.nx,grid.dim) !=
                       sub2ind(cindex,grid.nx,grid.dim))
   {
      SparseElt2 *temp;
      temp = new SparseElt2;
      (*temp).cindex = new int[grid.dim];
      for (r = 0; r < grid.dim; r++) 
         (*temp).cindex[r] = cindex[r];
      (*temp).val = value;
      (*temp).next = current;
      if (prev != NULL)
         (*prev).next = temp;
      else
         setvalarray(A,rindex,temp);
   }
   else
      (*prev).val += value;
}

void sparseorder(int *cindex, SparseElt2* &A, double value, GridData &grid)
{
   int r;
   SparseElt2 *current, *prev;

   for (prev = NULL,current = A;
        current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) <=
                           sub2ind(cindex,grid.nx,grid.dim); 
        prev = current,current = (*current).next);
   if (prev == NULL || sub2ind((*prev).cindex,grid.nx,grid.dim) !=
                       sub2ind(cindex,grid.nx,grid.dim))
   {
      SparseElt2 *temp;
      temp = new SparseElt2;
      (*temp).cindex = new int[grid.dim];
      for (r = 0; r < grid.dim; r++) 
         (*temp).cindex[r] = cindex[r];
      (*temp).val = value;
      (*temp).next = current;
      if (prev != NULL)
         (*prev).next = temp;
      else
         A = temp;
   }
   else
      (*prev).val += value;
}

void multsparserow2(SparseElt2 ****A, int *index, double value)
{
   SparseElt2 *current;

   for (current = evalarray(A,index); current != NULL; current = (*current).next)
      (*current).val *= value;
}

void multsparserow2(SparseElt2 *A, double value)
{
   SparseElt2 *current;

   for (current = A; current != NULL; current = (*current).next)
      (*current).val *= value;
}

void clearsparse(SparseElt2* &A)
{
   SparseElt2 *current;

   if (A != NULL)
   {
      for (current = (*A).next; current != NULL; current = (*current).next)
      {
         delete [] (*A).cindex;
         delete A;
         A = current;
      }
      delete [] (*A).cindex;
      delete A;
      A = NULL;
   }
}

void clearsparse(SparseElt2**** &A, GridData &grid)
{
   SparseElt2 *current1, *current2;
   int i, tindex[grid.dim];

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      current1 = evalarray(A,tindex);
      if (current1 != NULL)
      {
         for (current2 = (*current1).next; current2 != NULL; 
              current2 = (*current2).next)
         {
            delete [] (*current1).cindex;
            delete current1;
            current1 = current2;
         }
         delete [] (*current1).cindex;
         delete current1;
         setvalarray(A,tindex,NULL);
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}

void outputsparserow2(SparseElt2 *R, GridData &grid)
{
   SparseElt2 *current;

   if (R != NULL)
      for (current = R; current != NULL; current = (*current).next)
         cout << (*current).cindex[0] << " " << (*current).cindex[1] << " "
              << (*current).cindex[2] << " " << (*current).val << endl;
   else
      cout << "empty" << endl;
}

void outputsparsesmall(ofstream& outfile, SparseElt2 ****A, double ***S, PBData &pb, 
                       GridData &grid)
{
   int i, m, n, count = 0, tindex[grid.dim], rindex[grid.dim];
   SparseElt2 *current;
   double ehere;

   // outfile.precision(16);
   
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      count = 2*grid.dim+1;
      if (evalarray(S,tindex) < 0.0)
         ehere = pb.epsilonm;
      else
         ehere = pb.epsilonp;

      if (evalarray(A,tindex) == NULL)
      {
         for (m = 0; m < grid.dim; m++)
            rindex[m] = tindex[m];
         // outfile << scientific << sub2ind(tindex,grid.nx,grid.dim)+1 << " " 
         //         << sub2ind(tindex,grid.nx,grid.dim)+1 << " " 
         //         << 2.0*grid.dim*ehere/(grid.dx[0]*grid.dx[0]) << endl;
//         cout << sub2ind(tindex,grid.nx,grid.dim)+1 << " " 
//                 << sub2ind(tindex,grid.nx,grid.dim)+1 << " " 
//                 << 2.0*grid.dim*ehere/(grid.dx[m]*grid.dx[m]) << endl;
         for (m = 0; m < grid.dim; m++)
         {
            for (n = -1; n <= 1; n += 2)
            {
               rindex[m] = tindex[m]+n;
               if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
               {
                  // outfile << scientific << sub2ind(tindex,grid.nx,grid.dim)+1 << " " 
                  //         << sub2ind(rindex,grid.nx,grid.dim)+1 << " " 
                  //         << -ehere/(grid.dx[m]*grid.dx[m]) << endl;
//                  cout << sub2ind(tindex,grid.nx,grid.dim)+1 << " " 
//                          << sub2ind(rindex,grid.nx,grid.dim)+1 << " " 
//                          << -ehere/(grid.dx[m]*grid.dx[m]) << endl;
               }
            }
            rindex[m] = tindex[m];
         }
      }
      else
      {
//         count = 0;
//         for (current = evalarray(A,tindex); current != NULL; 
//              current = (*current).next)
//            count++;
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
         {
            // outfile << scientific << sub2ind(tindex,grid.nx,grid.dim)+1 << " "
            //         << sub2ind((*current).cindex,grid.nx,grid.dim)+1 << " "
            //         << (*current).val << endl;
//            if (count > 1)
//               cout << sub2ind(tindex,grid.nx,grid.dim)+1 << " "
//                       << sub2ind((*current).cindex,grid.nx,grid.dim)+1 << " "
//                       << (*current).val << endl;
         }
      }
//      if (count > 2*grid.dim+1)
//         getchar();
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}


// used to store matrix in iterative matrix solver
// fdstatus is array of size Nfd, initialized to be 0 
// fdstatus == 0 means no memory allocated; 1 means memory allocated and temporary;
// 2 means memory already in use
double ***setfourd(double ****fourd, char *fdstatus, int Nfd, GridData &grid)
{


   int i;

// looking for temporary, find fisrt i such that fdstatus[i]!=2, allocate matrix
   for (i = 0; i < Nfd && fdstatus[i] == 2; i++);
   if (i < Nfd)
   {
// if one is found, allocate memory if needed and fix
      if (fdstatus[i] == 0)
         fourd[i] = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
      fdstatus[i] = 2;

      return fourd[i];
   }
   else
   {
      cout << "Need more storage for fourd" << endl;
      exit(1);
   }
}

void removefourd(double ***r, double ****fourd, char *fdstatus, int Nfd)
{
   int i;

   for (i = 0; i < Nfd; i++)
      if (r != NULL && r == fourd[i])
         fdstatus[i] = 1;
}

void clearfourd(double**** &fourd, char* &fdstatus, int Nfd, GridData &grid)
{
   int i;

   for (i = 0; i < Nfd; i++)
      if (fourd[i] != NULL)
        free_matrix(fourd[i],grid.nx[0],grid.nx[1],grid.nx[2]);
   
   delete [] fourd;
   delete [] fdstatus;
}



void sparseorder(int row, int col, SparseElt** &A, double value)
{
   int r;
   SparseElt *current, *prev;

   for (prev = NULL,current = A[row]; current != NULL && (*current).cindex <= col;
        prev = current,current = (*current).next);
   if (prev == NULL || (*prev).cindex != col)
   {
      SparseElt *temp;
      temp = new SparseElt;
      (*temp).cindex = col;
      (*temp).val = value;
      (*temp).next = current;
      if (prev != NULL)
         (*prev).next = temp;
      else
         A[row] = temp;
   }
   else
      (*prev).val += value;
}

// used in AMGsmall
// A is the left most node for a row corresponding to A[tindex]
// roc is column index,
void sparseorder(int roc, SparseElt* &A, double value)
{
   int r;
   SparseElt *current, *prev;
   // start at A, stop if current is NULL or current.cindex > roc
   // prev = NULL and current = NULL : first time create diagonal node

   for (prev = NULL,current = A; current != NULL && (*current).cindex <= roc;
        prev = current,current = (*current).next);
   if (prev == NULL || (*prev).cindex != roc)
   {
      SparseElt *temp;
      temp = new SparseElt;
      (*temp).cindex = roc;
      (*temp).val = value;
      (*temp).next = current;
      if (prev != NULL)
         (*prev).next = temp; // prev!=Null, prev.cindex!= roc, current = NULL or current.cindex > roc : append 
      else
         A = temp; // prev==NULL, current = NULL or current.cindex > roc : prepend
   }
   else
      (*prev).val += value; // prev !=NULL and prev.cindex != roc: increment value
}

void removeelt(SparseElt* &current, SparseElt* &head)
{
   SparseElt *current2;

   if (current != NULL)
   {
      if (current != head)
      {
         for (current2 = head; (*current2).next != NULL && (*current2).next != current; 
              current2 = (*current2).next);
         if ((*current2).next == current)
            (*current2).next = (*current).next;
      }
      else
         head = (*current).next;

      delete current;
      current = NULL;
   }
}



void leftmultmxsmall(double ***y, SparseElt2 ****A, double ***x, GridData &grid,
                     double ***S, PBData &pb)
{
   int i, m, n, tindex[grid.dim], rindex[grid.dim];
   double ehere, temp;
   SparseElt2 *current;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (m = 0; m < grid.dim; m++)
         rindex[m] = tindex[m];
      if (evalarray(S,tindex) < 0.0)
         ehere = pb.epsilonm;
      else
         ehere = pb.epsilonp;

      temp = 0.0;
      if (evalarray(A,tindex) == NULL)
         for (m = 0; m < grid.dim; m++)
         {
            temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            for (n = -1; n <= 1; n += 2)
            {
               rindex[m] = tindex[m]+n;
               if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                  temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            }
            rindex[m] = tindex[m];
         }
      else
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
            temp += (*current).val*evalarray(x,(*current).cindex);

      setvalarray(y,tindex,temp);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}



void leftmultILUinv(double ***y, SparseElt2**** &M, double ***x, GridData &grid)
{
   int i, tindex[grid.dim];
   double d;
   SparseElt2 *current;
   int N = 1;
   for (i = 0; i < grid.dim; i++)
      N *= grid.nx[i]+1;
/*
   clock_t cstart1, cend1;
   clock_t cstart2, cend2;
   clock_t cstart3, cend3;
   clock_t cstart4, cend4;
   clock_t cstart5, cend5;
   clock_t cstart6, cend6;
   clock_t cstart7, cend7;
   clock_t cstart8, cend8;
   double temp;
   
   cstart1 = clock ();
   for (i = 0; i < N; i++)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
      setvalarray(y,tindex,0.0);
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
         setvalarray(y,tindex,evalarray(y,tindex)+
                              (*current).val*evalarray(x,(*current).cindex));
   }
   cend1 = clock ();
   cout << "mx-vec mult time = " << (double) (cend1-cstart1)/CLOCKS_PER_SEC << endl;
   cstart3 = clock ();
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(y,tindex,0.0);
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
         setvalarray(y,tindex,evalarray(y,tindex)+
                              (*current).val*evalarray(x,(*current).cindex));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend3 = clock ();
   cout << "mx-vec mult time = " << (double) (cend3-cstart3)/CLOCKS_PER_SEC << endl;
   cstart4 = clock ();
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = 0.0;
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
         temp += (*current).val*evalarray(x,(*current).cindex);
      setvalarray(y,tindex,temp);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend4 = clock ();
   cout << "mx-vec mult time = " << (double) (cend4-cstart4)/CLOCKS_PER_SEC << endl;
   cstart6 = clock ();
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(x,tindex);
      for (current = evalarray(M,tindex); 
           current != NULL && 
           ((*current).cindex[0] < tindex[0] || 
            ((*current).cindex[0] == tindex[0] && (*current).cindex[1] < tindex[1]) || 
            ((*current).cindex[0] == tindex[0] && (*current).cindex[1] == tindex[1] &&
             (*current).cindex[2] < tindex[2]));
           current = (*current).next)
         temp -= (*current).val*evalarray(y,(*current).cindex);
      setvalarray(y,tindex,temp);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend6 = clock ();
   cout << "forward solve time = " << (double) (cend6-cstart6)/CLOCKS_PER_SEC << endl;
   cstart7 = clock ();
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(x,tindex);
      for (current = evalarray(M,tindex); 
           current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) < 
                              sub2ind(tindex,grid.nx,grid.dim);
           current = (*current).next)
         temp -= (*current).val*evalarray(y,(*current).cindex);
      setvalarray(y,tindex,temp);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend7 = clock ();
   cout << "forward solve time = " << (double) (cend7-cstart7)/CLOCKS_PER_SEC << endl;
   int r;
   cstart8 = clock ();
   r = 0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(x,tindex);
      for (current = evalarray(M,tindex); 
           current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) < r;
           current = (*current).next)
         temp -= (*current).val*evalarray(y,(*current).cindex);
      setvalarray(y,tindex,temp);
      r++;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend8 = clock ();
   cout << "forward solve time = " << (double) (cend8-cstart8)/CLOCKS_PER_SEC << endl;
*/
//   cstart2 = clock ();
   for (i = 0; i < N; i++)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
      setvalarray(y,tindex,evalarray(x,tindex));
      for (current = evalarray(M,tindex); 
           current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) < i;
           current = (*current).next)
         setvalarray(y,tindex,evalarray(y,tindex)-
                              (*current).val*evalarray(y,(*current).cindex));
   }
//   cend2 = clock ();
//   cout << "forward solve time = " << (double) (cend2-cstart2)/CLOCKS_PER_SEC << endl;
//   cstart5 = clock ();
   for (i = N-1; i >= 0; i--)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
      for (current = evalarray(M,tindex); 
           current != NULL && sub2ind((*current).cindex,grid.nx,grid.dim) < i;
           current = (*current).next);
      d = (*current).val;
      for (current = (*current).next; current != NULL; current = (*current).next)
         setvalarray(y,tindex,evalarray(y,tindex)-
                              (*current).val*evalarray(y,(*current).cindex));
      setvalarray(y,tindex,evalarray(y,tindex)/d);
   }
//   cend5 = clock ();
//   cout << "backward solve time = " << (double) (cend5-cstart5)/CLOCKS_PER_SEC << endl;
//   exit(1);
}
