#include "storage.h"
#include "matrix.h"

#include "stdio.h"

extern int globsmall;
using namespace std;

void newstorage(StorageStruct &Dusmall)
{
   Dusmall.N = 4;
   Dusmall.info = new int[Dusmall.N];
   Dusmall.head = NULL;
}
// get initial storage size for Dusmall:  4 x number of points near interface
int getstoragesize(double ***S, GridData &grid)
{
   int count = 0, i, r, s, tindex[grid.dim], rindex[grid.dim];

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = tindex[r]+s;
            if (rindex[r] >= 0 && rindex[r] <= grid.nx[r] && 
                (evalarray(S,tindex) <= 0.0)+(evalarray(S,rindex) <= 0.0) == 1)
               count++;//if tindex has neighbor(among 6) across interface
         }
         rindex[r] = tindex[r];
      }
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

// 3 for Du and 1 for uint
   return 4*count;
}

int findinstorage(int *info, StorageStruct *Dusmall, int smallsize)
{
   int r;

   if (smallsize == 1)
   {
      for (r = 0; r < Dusmall[0].N && info[r] == Dusmall[0].info[r]; r++);
      if (r >= Dusmall[0].N)
         return 0;
      else
      {
         cout << "Error: not in storage " << info[0] << " " << info[1] << " " 
              << info[2] << " " << info[3] <<  endl;
         exit(1);
      }
   }
   else
   {
      int loc, halfsize = (int) (smallsize-1)/2;

      for (r = 0; r < Dusmall[halfsize].N && info[r] == Dusmall[halfsize].info[r]; r++);
      if (r == Dusmall[halfsize].N || info[r] <= Dusmall[halfsize].info[r])
         loc = findinstorage(info,Dusmall,halfsize+1); 
      else
      {
         loc = findinstorage(info,&(Dusmall[halfsize+1]),smallsize-(halfsize+1)); 
         loc += halfsize+1;
      }
/*
      for (r = 0; r < Dusmall[loc].N && info[r] == Dusmall[loc].info[r]; r++);
      if (r >= Dusmall[loc].N);
      else
      {
         cout << "Error: not in storage " << info[0] << " " << info[1] << " " 
              << info[2] << " " << info[3] <<  endl;
         exit(1);
      }
*/

      return loc;
   }
}

void evalfromstorage(double &uint, double *Du, int *index, int rstar, int sstar,
                     int mid, StorageStruct *Dusmall, int smallsize, double ***u, 
                     double ***S, PBData &pb, GridData &grid)
{
// just like getinterfaceDu
   int i, r, t, loc, N = 2*mid; 
   int info[Dusmall[0].N], Narray[grid.dim], sindex[grid.dim], tindex[grid.dim];
   double value;
   SparseElt *current;

   for (r = 0; r < grid.dim; r++)
      Narray[r] = N;

   info[0] = sub2ind(index,grid.nx,grid.dim);
   info[1] = rstar;
   info[2] = sstar;

   for (t = -1; t < grid.dim; t++)
   {
      info[3] = t;
      value = 0.0;

      loc = findinstorage(info,Dusmall,smallsize);
//      cout << info[0] << " " << info[1] << " " << info[2] << " " << info[3] << endl;
//      cout << "loc = " << loc << endl;

      for (current = Dusmall[loc].head; current != NULL; current = (*current).next)
      {
         if ((*current).cindex >= 0)
         {
            if (globsmall == 1)
            {
               ind2sub(sindex,(*current).cindex,Narray,grid.dim);
               for (r = 0; r < grid.dim; r++)
                  tindex[r] = index[r]+sindex[r]-mid;
            }
            else if (globsmall == 2)
               ind2sub(tindex,(*current).cindex,grid.nx,grid.dim);
            value += (*current).val*evalarray(u,tindex);
         }
         else
            value += (*current).val;
      }

      if (t < 0)
         uint = value;
      else
         Du[t] = value;
   }
}
// no parameter for mid, used in cim
void evalfromstorage(double &uint, double *Du, int *index, int rstar, int sstar,
                     StorageStruct *Dusmall, int smallsize, double ***u, double ***S, 
                     PBData &pb, GridData &grid)
{
// just like getinterfaceDu
   int i, r, t, loc, N, mid;
   int info[Dusmall[0].N], Narray[grid.dim], sindex[grid.dim], tindex[grid.dim];
   double value;
   SparseElt *current;

   info[0] = sub2ind(index,grid.nx,grid.dim);
   info[1] = rstar;
   info[2] = sstar;

   for (t = -1; t < grid.dim; t++)
   {
      info[3] = t;
      value = 0.0;

      loc = findinstorage(info,Dusmall,smallsize);

      for (current = Dusmall[loc].head; current != NULL; current = (*current).next)
      {
         if ((*current).cindex >= 0)
         {
            if (globsmall == 1)
            {
               mid = Dusmall[loc].mid;
               N = 2*mid;
               for (r = 0; r < grid.dim; r++)
                  Narray[r] = N;
               ind2sub(sindex,(*current).cindex,Narray,grid.dim);
               for (r = 0; r < grid.dim; r++)
                  tindex[r] = index[r]+sindex[r]-mid;
            }
            else if (globsmall == 2)
               ind2sub(tindex,(*current).cindex,grid.nx,grid.dim);
            value += (*current).val*evalarray(u,tindex);
         }
         else
            value += (*current).val;
      }

      if (t < 0)
         uint = value;
      else
         Du[t] = value;
   }
}

void clearstorage(StorageStruct* &Dusmall, int &smallsize)
{
   int i;
   SparseElt *current, *current2;

   for (i = 0; i < smallsize; i++)
   {
      for (current = Dusmall[i].head, current2 = current;
           current != NULL;
           current = (*current).next, delete current2,
           current2 = current)
         delete [] Dusmall[i].info;
   }
   delete [] Dusmall;

   smallsize = 0;
}
