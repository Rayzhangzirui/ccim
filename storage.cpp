#include "storage.h"
#include "matrix.h"
#include "stdio.h"

#include "interface.h"
#include "global.h"

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


// look up info in array of Dusmall using binary search
// info  = [ind, r, sk, t] all are ordered
// ind is ordered, 0 to N^3, r from 0 to 2, sk -1 1, t 0 to 2 
// Dusmall[s].N = 4 
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

      // if r == 4, then found match
      // 
      if (r == Dusmall[halfsize].N || info[r] <= Dusmall[halfsize].info[r])
         loc = findinstorage(info, Dusmall, halfsize+1); 
      else
      {
         loc = findinstorage(info,&(Dusmall[halfsize+1]),smallsize-(halfsize+1)); 
         loc += halfsize+1;
      }
      return loc;
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
   SparseElt *current, *current2;
   for (int i = 0; i < smallsize; i++)
   {
      for (current = Dusmall[i].head, current2 = current; current != NULL; 
           current = (*current).next, delete current2, current2 = current);
      delete [] Dusmall[i].info;
   }
   delete [] Dusmall;
   smallsize = 0;
}



// check du with StorageStruct at interface
void checkDuStorage(double ***u, StorageStruct *Dusmall, int smallsize, double ***S, PBData &pb, GridData &grid)
{
  cout<< "[checkcim345Du Storage]" <<endl;
   int i, r, s, t, m, mid = 2;
   int tindex[grid.dim], rindex[grid.dim];
   double uint, Du[grid.dim];
   double alpha, thesign, tangent[grid.dim], normal[grid.dim];
   double tmperr;
   
   int zindex[grid.dim][grid.dim];

   double theerr[grid.dim];
   double inside_max_err[grid.dim];
   double outside_max_err[grid.dim];

   for (t = 0; t < grid.dim; t++){
      theerr[t] = 0.0;
      inside_max_err[t] = 0.0;
      outside_max_err[t] = 0.0;
   }

   double uint2, Du2[grid.dim];

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (evalarray(S,tindex) < 0.0)
         thesign = -1.0;
      else
         thesign = 1.0;
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = tindex[r]+s;
            if (rindex[r] >= 0 && rindex[r] <= grid.nx[r] && 
                (evalarray(S,tindex) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)
            {
//               getinterfaceDu(uint2,Du2,tindex,r,s,u,S,pb,grid);
//               evalfromstorage(uint,Du,tindex,r,s,mid,Dusmall,smallsize,u,S,pb,grid);
               evalfromstorage(uint,Du,tindex,r,s,Dusmall,smallsize,u,S,pb,grid);

               getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
               for (t = 0; t < grid.dim; t++)
               {
                  tmperr = fabs(Du[t]-getDu(tindex,t,r,s,alpha,thesign,grid));
                  if (tmperr > theerr[t])
                  {
                     theerr[t] = tmperr;
                     for (m = 0; m < grid.dim; m++)
                        zindex[t][m] = tindex[m];
                  }

                  if (tmperr > inside_max_err[t] && thesign < 0.0){
                     inside_max_err[t] = tmperr;
                  }
                     
                  if (tmperr > outside_max_err[t] && thesign + 0.0){
                     outside_max_err[t] = tmperr;
                  }
               }

               // if (globwriteerr){
               //    outfile_Duerr<<tindex[0]<<","<<tindex[1]<<","<<tindex[2]<<","<<r<<","<<s<<",";
               //    outfile_Duerr <<setprecision(12)
               //    <<Du[0]-getDu(tindex,0,r,s,alpha,thesign,grid)<<","
               //    <<Du[1]-getDu(tindex,1,r,s,alpha,thesign,grid)<<","
               //    <<Du[2]-getDu(tindex,2,r,s,alpha,thesign,grid)
               //    <<endl;
               //  }

            }
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
   
   cout << "Error in cim345 derivative is " << theerr[0] << " " << theerr[1] << " " 
        << theerr[2] << endl;
   cout << "   at";
   for (t = 0; t < grid.dim; t++)
   {
      for (m = 0; m < grid.dim; m++)
         cout << " " << zindex[t][m];
      // cout << " (" << getstatus5(S,zindex[t],grid) << ")";
      if (t < grid.dim-1)
         cout << ", ";
      else
         cout << "." << endl;
   }

   cout << "Inside max err is " << inside_max_err[0] << " " << inside_max_err[1] << " " << inside_max_err[2] << endl;
   cout << "Outside max err is " << outside_max_err[0] << " " << outside_max_err[1] << " "<< outside_max_err[2] << endl;

}

