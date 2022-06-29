#include "march.h"
#include "matrix.h"
#include "global.h"
#include "getvn.h"

#include <iostream>
#include <cmath>


using namespace std;

// initialize marching structure, S is surface
void init_march(MarchStruct &march, double*** S, PBData &pb, GridData&grid){
   march.dist = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   march.status = cmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   march.extend = new double ***[1];

   march.extend[0] = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   march.nval = 1;

   march.dorig = S; // surface

   march.tol = 1.0e-12;
   (march.heap).head = NULL;
   (march.heap).tail = NULL;
   (march.heap).loc = heapmatrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   (march.heap).dim = grid.dim;
   march.pb = pb;

   march.smallsize = 0;
}

// fast marching
void fmarch(MarchStruct &march, double phitube, GridData &grid)
{
   int i, r, s;
   int tindex[grid.dim], rindex[grid.dim];
   char notstop;
   clock_t cstart, cend;

   // === go over every grid point, set status as FAR(1)
   for (i = 0; i < grid.dim; i++) 
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(march.status,tindex,1);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   // 
   cstart = clock ();
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      // if status is not ACCEPT(3)
      if (evalarray(march.status,tindex) != 0 && evalarray(march.status,tindex) != 3)
      {
         if (evalarray(march.dorig,tindex) != 0.0)
            notstop = 1; //if tindex not exactly at interface
         else
            notstop = 0;//if tindex exactly at interface, should stop

         // for points not exactly on interface, if any nbr is on diff side, mark as stop
         for (r = 0; r < grid.dim; r++)
            rindex[r] = tindex[r];
         for (r = 0; r < grid.dim && notstop; r++)
         {
            for (s = -1; s <= 1 && notstop; s += 2)
            {
               rindex[r] = tindex[r]+s;
               if (rindex[r] >= 0 && rindex[r] <= grid.nx[r] && 
                   evalarray(march.status,rindex) != 0 && 
                   (evalarray(march.dorig,tindex) <= 0.0)+
                   (evalarray(march.dorig,rindex) <= 0.0) == 1)
                  notstop = 0;
            }
            rindex[r] = tindex[r];
         }

         // for points near interface
         if (!notstop) 
         {
            fmarchstep(march,tindex,grid); //get temp distance and extension
            addtoheap(march.heap,tindex,march.dist);//add to heap
            setvalarray(march.status,tindex,2); //set status to CLOSE
         }
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cend = clock ();
   cout << "fmarch initialization = " << (double) (cend-cstart)/CLOCKS_PER_SEC 
        << endl;


   while ((march.heap).head != NULL && (*((march.heap).head)).index[0] >= 0)
      if (phitube < 0.0 ||
          fabs(evalarray(march.dist,(*((march.heap).head)).index)) <= phitube)
      {
         setvalarray(march.status,(*((march.heap).head)).index,3); // set status ACCEPT(3) head of heap, smallest distance
         for (i = 0; i < grid.dim; i++)
            tindex[i] = (*((march.heap).head)).index[i]; // get index of head
         fixheapeltdelete(march.heap,(march.heap).head,march.dist); // remove head from heap

         // go through nbhr of tindex and modify distance
         for (r = 0; r < grid.dim; r++) 
            rindex[r] = tindex[r];
         for (r = 0; r < grid.dim; r++)
         {
            for (s = -1; s <= 1; s += 2)
            {
               rindex[r] = tindex[r]+s;
               if (rindex[r] >= 0 && rindex[r] <= grid.nx[r] && 
                   evalarray(march.status,rindex) != 0 && 
                   evalarray(march.status,rindex) != 3)
               {
                   fmarchstep(march,rindex,grid);// calculate distance and extension of rindex
                   if (evalarray(march.status,rindex) == 1) // if status is FAR(1)
                   {
                      addtoheap(march.heap,rindex,march.dist);// add to heap
                      setvalarray(march.status,rindex,2); // set status CLOSE(2)
                   }
                   else // if status is CLOSE(2), modify hep
                      fixheapeltreplace(march.heap,rindex,march.dist);
               }
            }
            rindex[r] = tindex[r];
         }
      }
      else
         fixheapeltdelete(march.heap,(march.heap).head,march.dist);


   // if fast marhing to bandwidth phitube, set phi>phitube to (+/-)phitube
   if (phitube >= 0.0) 
   {
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         if (evalarray(march.status,tindex) != 3)
         {
            if (evalarray(march.dorig,tindex) <= 0.0)
               setvalarray(march.dist,tindex,-phitube);
            else
               setvalarray(march.dist,tindex,phitube);
         }

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
}


// update distance and extension at index with nbr information
void fmarchstep(MarchStruct &march, int *index, GridData &grid)
{
   int r;
   double dvalue[grid.dim];
   double **evalue = matrix(march.nval-1,grid.dim-1); //1 by 3 matrix
   double thedx[grid.dim];
   char yes[grid.dim];

   if (evalarray(march.status,index) != 0 && evalarray(march.status,index) != 3)
   {
      if (evalarray(march.dorig,index) != 0.0)
      {
         fmarchdirection(yes,dvalue,evalue,thedx,march,index,grid);
         while (!fmarchdistance(march,yes,dvalue,evalue,thedx,index,grid));
      }   
      else
      {
         cout << "This fmarch option should not happen" << endl;
         exit(1);
         setvalarray(march.dist,index,evalarray(march.dorig,index));
         for (r = 0; r < march.nval; r++)
// CHANGE HERE
// insert your program for calculating normal velocity
            setvalarray(march.extend[r],index,getvn((march.pb).psi,(march.pb).LJ,
                        march.dorig,index,0,march.pb,grid));
// END CHANGE
      }
   }

   free_matrix(evalue,march.nval-1,grid.dim-1);
}

// get direction to look for update, if yes[r] = 1, use info comming in dim r, dvalue[r] = phi value in dim r, 
// evalue[0][r] = vn value in dim r, thedx[r] = dx in dim r
void fmarchdirection(char *yes, double *dvalue, double **evalue, double *thedx,
                     MarchStruct &march, int *index, GridData &grid)
{
   char yesper[2]; // yesper[0] = 1 if -1 direction change sign, 0 if -1 direction same sign. yesper[1] for +1 direction
   double phiper[2], valper[march.nval][2], dxper[2], x[grid.dim];
   int rindex[grid.dim];
   int j, r, s;

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
   {
      for (s = 0; s <= 1; s++)
         yesper[s] = 0;
      for (s = -1; s <= 1; s += 2)
      {
         rindex[r] = index[r]+s;
         if (rindex[r] >= 0 && rindex[r] <= grid.nx[r])
         {
            if ((evalarray(march.dorig,index) <= 0.0)+
                (evalarray(march.dorig,rindex) <= 0.0) != 1 &&
                evalarray(march.dorig,index) != 0.0) //if index away from interface
            {
               if (evalarray(march.status,rindex) == 3) //if rindex already ACCEPT(3)
               {
                  yesper[(s+1)/2] = 1;
                  phiper[(s+1)/2] = evalarray(march.dist,rindex);
                  dxper[(s+1)/2] = grid.dx[r];
                  for (j = 0; j < march.nval; j++)
                     valper[j][(s+1)/2] = evalarray(march.extend[j],rindex);
               }
            }
            else//if index near interface
            {
               yesper[(s+1)/2] = 1;
               phiper[(s+1)/2] = 0.0;
               if (evalarray(march.dorig,index) != 0.0)
                  dxper[(s+1)/2] = -evalarray(march.dorig,index)*grid.dx[r]/
                                    (evalarray(march.dorig,rindex)-
                                     evalarray(march.dorig,index)); // distance to interface in dim r
               else
                  dxper[(s+1)/2] = 0.0;
               for (j = 0; j < march.nval; j++)
               {
// CHANGE HERE
// insert your program for calculating normal velocity
                  if (!globheat && !globexactvn)
                     if (s == 1)
                        valper[j][(s+1)/2] = getvn((march.pb).psi,(march.pb).LJ,
                                                   march.dorig,index,r,march.pb,grid);
                     else
                        valper[j][(s+1)/2] = getvn((march.pb).psi,(march.pb).LJ,
                                                   march.dorig,rindex,r,march.pb,grid);
                  else if (globheat && !globexactvn)
                     if (!globsmall)
                        valper[j][(s+1)/2] = getheatvn((march.pb).psi,march.dorig,index,
                                                       r,s,march.pb,grid);
                     else
                        valper[j][(s+1)/2] = getheatvn((march.pb).psi,march.dorig,index,
                                                       r,s,march.Dusmall,
                                                       march.smallsize,march.pb,grid);
                  else if (globexactvn)
                     valper[j][(s+1)/2] = getexactvn(index,r,s,march.dorig,grid);
/*
                  if (index[0] == grid.nx[0]/2 && index[1] == grid.nx[2]/2 && 
                      index[2] >= grid.nx[0]/2 && r == 2 && s == 1)
                  {
                     double alpha, tangent[grid.dim], normal[grid.dim], Dup[grid.dim],
                            Dum[grid.dim];
                     int t;
                     getinterfaceinfo(alpha,tangent,normal,march.dorig,index,r,s,grid);
                     for (t = 0; t < grid.dim; t++)
                     {
                        Dup[t] = getDu(index,t,r,s,alpha,1.0,grid);
                        Dum[t] = getDu(index,t,r,s,alpha,-1.0,grid);
                     }
                     mu2d.precision(16);
                     mu2d << scientific << grid.t << " " << valper[j][(s+1)/2] << " " 
                          << (1.0-(march.pb).epsilonp/(march.pb).epsilonm)*grid.radius/
                             sqrt(grid.radius*grid.radius+
                                  fabs(1.0-(march.pb).epsilonp/
                                           (march.pb).epsilonm)) << " "
                          << formvn(Dup,Dum,normal,grid) << " "
                          << alpha << " "
                          << grid.a[2]+(index[2]+alpha)*grid.dx[2] << " "
                          << grid.radius << endl;
                  }
//                  mu2d << valper[j][(s+1)/2] << endl;
*/
//END CHANGE
               }
            }
         }
      }
      yes[r] = 0;
      dvalue[r] = 0.0;
      for (j = 0; j < march.nval; j++)
         evalue[j][r] = 0.0;
      thedx[r] = 0.0;
      for (s = 0; s <= 1; s++)
         if ((yesper[s] && !yesper[(s+1)%2]) || // s(=-1/+1) side ACCEPT, the other side not ACCEPT OR
             (yesper[s] && yesper[(s+1)%2] &&  //both side ACCEPT AND
// changed below line from <= to < since former is wrong
              (dxper[s] < dxper[(s+1)%2] || //dx on s side is smaller
               (dxper[s] == dxper[(s+1)%2] && // dx is the same on both side,  
                fabs(phiper[s]) <= fabs(phiper[(s+1)%2]))))) // phi on s side is smaller
         {
            dvalue[r] += phiper[s];
            for (j = 0; j < march.nval; j++)
               evalue[j][r] += valper[j][s];
            thedx[r] = dxper[s];
            if (yes[r])
            {
               dvalue[r] /= 2.0;
               for (j = 0; j < march.nval; j++)
                  evalue[j][r] /= 2.0;
            }
            else
               yes[r] = 1;
         }
      rindex[r] = index[r];
   }
}
// get the distance and extension using nbhr info, dvalue for phi, evalue for vn, thedx for distance
char fmarchdistance(MarchStruct &march, char *yes, double *dvalue, double **evalue,
                    double *thedx, int *index, GridData &grid)
{
   double a = 0.0, b = 0.0, c = -1.0, maxvalue, numerator[march.nval], denominator;
   int i, imax = 0, r;
   char gotit = 0;

   maxvalue = 0.0;
   for (i = 0; i < grid.dim; i++)
      if (yes[i] && fabs(dvalue[i]) > maxvalue)
      {
         maxvalue = fabs(dvalue[i]);
         imax = i;
      }
   for (i = 0; i < grid.dim; i++)
      if (yes[i] && fabs(thedx[i]) < march.tol)
         thedx[i] = march.tol;
   for (i = 0; i < grid.dim && (!(yes[i]) || fabs(thedx[i]) <= march.tol); i++);

   if (i < grid.dim)
   {
      for (i = 0; i < grid.dim; i++)
         if (yes[i])
         {
            a += 1.0/(thedx[i]*thedx[i]);
            b += -2.0*dvalue[i]/(thedx[i]*thedx[i]);
            c += dvalue[i]*dvalue[i]/(thedx[i]*thedx[i]);
         }
      if (b*b-4.0*a*c >= 0.0)
      {
         if (evalarray(march.dorig,index) < 0.0)
            setvalarray(march.dist,index,(-b-sqrt(b*b-4.0*a*c))/(2.0*a));
         else
            setvalarray(march.dist,index,(-b+sqrt(b*b-4.0*a*c))/(2.0*a));
         for (i = 0; i < grid.dim &&
                     (!(yes[i]) ||
                      fabs(evalarray(march.dist,index)) >= fabs(dvalue[i])); i++);
         if (i < grid.dim) // if exist i such that yes[i]=true but dist < dvalue[i], i.e. the newly calculated phi is smaller than nbr phi
            yes[imax] = 0;// set dim i with maximum phi, redo calculation
         else // new phi is valid, calculate vn by dot(grad(vn), grad(phi)) = 0
         {
            denominator = 0.0;
            for (i = 0; i < grid.dim; i++)
               if (yes[i])
                  denominator += (evalarray(march.dist,index)-dvalue[i])/
                                 (thedx[i]*thedx[i]);
            for (r = 0; r < march.nval; r++)
            {
               numerator[r] = 0.0;
               for (i = 0; i < grid.dim; i++)
                  if (yes[i])
                     numerator[r] += evalue[r][i]*(evalarray(march.dist,index)-
                                                   dvalue[i])/(thedx[i]*thedx[i]);
               setvalarray(march.extend[r],index,numerator[r]/denominator);
            }
            gotit = 1;
         }
      }
      else// no solution, redo calculation without info from largest phi
         yes[imax] = 0;
   }
   else // if all dimension has no coming info or distances<tol
   {
      gotit = 1;
      numerator[0] = 0.0;
      denominator = 0.0;
      for (i = 0; i < grid.dim; i++)
         if (yes[i])
         {
            numerator[0] += dvalue[i];
            denominator += 1.0;
         }
      setvalarray(march.dist,index,numerator[0]/denominator); // set dist as average of small distances (<tol)
      for (r = 0; r < march.nval; r++)
      {
         numerator[r] = 0.0;
         for (i = 0; i < grid.dim; i++)
            if (yes[i])
               numerator[r] += evalue[r][i];
         setvalarray(march.extend[r],index,numerator[r]/denominator); // set vn extension as avaerge 
      }
   }

   return gotit;
}

void addtoheap(HeapStruct &heap, int *index, double ***phi)
{
   HeapElt *parent, *child;
   int j, iparent, ichild, last;

   if (heap.tail != NULL)
      last = (*(heap.tail)).num+1;
   else
      last = 1;
   iparent = last>>1;
   ichild = last%2;
   if (iparent > 0)
      parent = heapgoto(heap,iparent);
   else
      parent = NULL;
   if (parent != NULL)
      if ((*parent).child[ichild] == NULL)
      {
         (*parent).child[ichild] = new HeapElt;
         child = (*parent).child[ichild];
         (*child).parent = parent;
         (*child).child[0] = NULL;
         (*child).child[1] = NULL;
      }
      else
         child = (*parent).child[ichild];
   else
      if (heap.head == NULL)
      {
         heap.head = new HeapElt;
         child = heap.head;
         (*child).parent = NULL;
         (*child).child[0] = NULL;
         (*child).child[1] = NULL;
      }
      else
         child = heap.head;

   (*child).num = last;
   for (j = 0; j < heap.dim; j++)
      (*child).index[j] = index[j];
   setvalarray(heap.loc,(*child).index,child);
   heap.tail = child;

   fixheapeltup(heap,child,phi);
}

void fixheapeltup(HeapStruct &heap, HeapElt *fix, double ***phi)
{
   HeapElt *parent, *child, temp;
   char change = 1;
   int j;

   child = fix;
   parent = (*child).parent;
   while (change && parent != NULL)
   {
      if (fabs(evalarray(phi,(*parent).index)) > fabs(evalarray(phi,(*child).index)))
      {
         for (j = 0; j < heap.dim; j++)
         {
            temp.index[j] = (*parent).index[j];
            (*parent).index[j] = (*child).index[j];
            (*child).index[j] = temp.index[j];
         }
         setvalarray(heap.loc,(*parent).index,parent);
         setvalarray(heap.loc,(*child).index,child);
         child = parent;
         parent = (*parent).parent;
      }
      else
         change = 0;
   }
}

HeapElt* fixheapeltempty(HeapStruct &heap, HeapElt *fix, double ***phi)
{
   HeapElt *current, *child[2];
   int j, r;

   setvalarray(heap.loc,(*fix).index,NULL);
   current = fix;
   child[0] = (*current).child[0];
   child[1] = (*current).child[1];
   while (child[0] != NULL && (*(child[0])).index[0] >= 0)
   {
      if (child[1] == NULL || (*(child[1])).index[0] < 0 ||
          fabs(evalarray(phi,(*(child[0])).index)) <=
          fabs(evalarray(phi,(*(child[1])).index)))
         r = 0;
      else
         r = 1;
      for (j = 0; j < heap.dim; j++)
         (*current).index[j] = (*(child[r])).index[j];
      setvalarray(heap.loc,(*current).index,current);
      current = child[r];
      child[0] = (*current).child[0];
      child[1] = (*current).child[1];
   }

   return current;
}

void fixheapeltdelete(HeapStruct &heap, HeapElt *del, double ***phi)
{
   HeapElt *current;
   int j;

   current = fixheapeltempty(heap,del,phi);
   if (current != heap.tail)
   {
      for (j = 0; j < heap.dim; j++)
         (*current).index[j] = (*(heap.tail)).index[j];
      setvalarray(heap.loc,(*current).index,current);
      for (j = 0; j < heap.dim; j++)
         (*(heap.tail)).index[j] = -1;
      heap.tail = heapgoto(heap,(*(heap.tail)).num-1);
      fixheapeltup(heap,current,phi);
   }
   else
   {
      for (j = 0; j < heap.dim; j++)
         (*(heap.tail)).index[j] = -1;
      heap.tail = heapgoto(heap,(*(heap.tail)).num-1);
   }
}

void fixheapeltreplace(HeapStruct &heap, int *index, double ***phi)
{
   HeapElt *replace, *current, temp;
   int j;

   replace = evalarray(heap.loc,index);

   for (j = 0; j < heap.dim; j++)
      temp.index[j] = (*replace).index[j];

   current = fixheapeltempty(heap,replace,phi);
   for (j = 0; j < heap.dim; j++)
      (*current).index[j] = temp.index[j];
   setvalarray(heap.loc,(*current).index,current);
   fixheapeltup(heap,current,phi);
}

HeapElt* heapgoto(HeapStruct &heap, int num)
{
   HeapElt *current;
   int n, k;

   if (num <= 0)
      return NULL;

   for (n = 1; n <= num; n*= 2);
   n /= 2;

   current = heap.head;
   k = num;
   while (n > 1 && current != NULL)
   {
      if (k >= n)
         k -= n;
      n /= 2;
      if (k >= n)
         current = (*current).child[1];
      else
         current = (*current).child[0];
   }

   return current;
}

HeapElt* evalarray(HeapElt ***A, int *index)
{
   return A[index[0]][index[1]];
}

HeapElt* evalarray(HeapElt ****A, int *index)
{
   return A[index[0]][index[1]][index[2]];
}

void setvalarray(HeapElt ***A, int *index, HeapElt *value)
{
   A[index[0]][index[1]] = value;
}

void setvalarray(HeapElt ****A, int *index, HeapElt *value)
{
   A[index[0]][index[1]][index[2]] = value;
}

HeapElt ***heapmatrix(int row, int col)
{
   int i, j;
   HeapElt ***m;
   int err = 0;

   m = new HeapElt**[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new HeapElt*[col+1];
      if (!m[i] && !err)
         err = 1;
      for (j = 0; j <= col; j++)
         m[i][j] = NULL;
   }
   if (err)
      cout << "matrix stuff didn`t work" << endl;
   return m;
}

HeapElt ****heapmatrix(int row, int col, int frb)
{
   int i, j, k;
   HeapElt ****m;
   int err = 0;

   m = new HeapElt***[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new HeapElt**[col+1];
      if (!m[i] && !err)
         err = 1;
      for (j = 0; j <= col; j++)
      {
         m[i][j] = new HeapElt*[frb+1];
         if (!m[i][j] && !err)
            err = 1;
         for (k = 0; k <= frb; k++)
            m[i][j][k] = NULL;
      }
   }
   if (err)
      cout << "matrix stuff didn`t work" << endl;
   return m;
}

char checkcompat(HeapStruct heap)
{
   int i, count = 0;
   HeapElt *current = heap.head, *temp;

   for (i = 1; current != NULL && (*current).index[0] >= 0; i++)
   {
      current = heapgoto(heap,i);
      if (current != NULL && (*current).index[0] >= 0)
      {
         count++;
         temp = evalarray(heap.loc,(*current).index);
         if (current != temp)
         {
            cout << "problem at heap number " << i << " with " 
                 << (*current).index[0] << " " << (*current).index[1] << " "
                 << (*current).index[2] << endl;
            getchar();
            return 0;
         }
      }
   }

   cout << "heap compatible " << count << endl;

   return 1;
}

void outputheap(HeapStruct heap, double ***value, double ***ovalue)
{
   int i, count = 0;
   HeapElt *current = heap.head;

   cout << "start heap write" << endl;
   for (i = 1; current != NULL && (*current).index[0] >= 0; i++)
   {
      current = heapgoto(heap,i);
      if (current != NULL && (*current).index[0] >= 0)
      {
         count++;
         cout << (*current).index[0] << " " << (*current).index[1] << " "
              << (*current).index[2] << " " << evalarray(value,(*current).index) 
              << " " << evalarray(ovalue,(*current).index) << endl;
      }
   }
   cout << "end heap write" << endl;
}
