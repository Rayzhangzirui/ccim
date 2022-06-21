#include "solver.h"
#include <iostream>
#include <cmath>
using namespace std;


void copyfromsmall(SparseElt2**** &B, SparseElt2**** &A, GridData &grid, double ***S, 
                   PBData &pb)
{
   int i, r, m, n, tindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   double ehere;

// copy A to B and order B
   clearsparse(B,grid);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(B,tindex,NULL);
      if (evalarray(A,tindex) == NULL)
      {
// create rows in B when they don't exist in A
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
         temp = new SparseElt2;
         (*temp).cindex = new int[grid.dim];
         for (r = 0; r < grid.dim; r++) 
            (*temp).cindex[r] = tindex[r];
         (*temp).val = 0.0;
         for (m = 0; m < grid.dim; m++)
            (*temp).val += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
         (*temp).next = NULL;
         setvalarray(B,tindex,temp);

         for (m = 0; m < grid.dim; m++)
            for (n = -1; n <= 1; n += 2)
            {
               temp = new SparseElt2;
               (*temp).cindex = new int[grid.dim];
               for (r = 0; r < grid.dim; r++) 
                  (*temp).cindex[r] = tindex[r];
               (*temp).cindex[m] = tindex[m]+n;
               if ((*temp).cindex[m] >= 0 && (*temp).cindex[m] <= grid.nx[m])
               {
                  (*temp).val = -ehere/(grid.dx[m]*grid.dx[m]);
                  for (prev = NULL,current2 = evalarray(B,tindex); 
                       current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                           sub2ind((*temp).cindex,grid.nx,grid.dim); 
                       prev = current2,current2 = (*current2).next);
                  (*temp).next = current2;
                  if (prev != NULL)
                     (*prev).next = temp;
                  else
                     setvalarray(B,tindex,temp);
               }
            }
      }
      else
         for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
         {
            temp = new SparseElt2;
            (*temp).cindex = new int[grid.dim];
            for (r = 0; r < grid.dim; r++) 
               (*temp).cindex[r] = (*current).cindex[r];
            (*temp).val = (*current).val;
            for (prev = NULL,current2 = evalarray(B,tindex); 
                 current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                     sub2ind((*temp).cindex,grid.nx,grid.dim); 
                 prev = current2,current2 = (*current2).next);
            (*temp).next = current2;
            if (prev != NULL)
               (*prev).next = temp;
            else
               setvalarray(B,tindex,temp);
         }   
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}

void copyfromsmall(SparseElt2**** &B, SparseElt2**** &A, GridData &grid, double ***a, 
                   double ***S, PBData &pb)
{
   int i, r, m, n, tindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   double ehere;

// copy A to B and order B
   clearsparse(B,grid);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(B,tindex,NULL);
      if (evalarray(A,tindex) == NULL)
      {
// create rows in B when they don't exist in A
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
         temp = new SparseElt2;
         (*temp).cindex = new int[grid.dim];
         for (r = 0; r < grid.dim; r++) 
            (*temp).cindex[r] = tindex[r];
         (*temp).val = evalarray(a,tindex);
         for (m = 0; m < grid.dim; m++)
            (*temp).val += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
         (*temp).next = NULL;
         setvalarray(B,tindex,temp);

         for (m = 0; m < grid.dim; m++)
            for (n = -1; n <= 1; n += 2)
            {
               temp = new SparseElt2;
               (*temp).cindex = new int[grid.dim];
               for (r = 0; r < grid.dim; r++) 
                  (*temp).cindex[r] = tindex[r];
               (*temp).cindex[m] = tindex[m]+n;
               if ((*temp).cindex[m] >= 0 && (*temp).cindex[m] <= grid.nx[m])
               {
                  (*temp).val = -ehere/(grid.dx[m]*grid.dx[m]);
                  for (prev = NULL,current2 = evalarray(B,tindex); 
                       current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                           sub2ind((*temp).cindex,grid.nx,grid.dim); 
                       prev = current2,current2 = (*current2).next);
                  (*temp).next = current2;
                  if (prev != NULL)
                     (*prev).next = temp;
                  else
                     setvalarray(B,tindex,temp);
               }
            }
      }
      else
         for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
         {
            temp = new SparseElt2;
            (*temp).cindex = new int[grid.dim];
            for (r = 0; r < grid.dim; r++) 
               (*temp).cindex[r] = (*current).cindex[r];
            (*temp).val = (*current).val;
            for (prev = NULL,current2 = evalarray(B,tindex); 
                 current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                     sub2ind((*temp).cindex,grid.nx,grid.dim); 
                 prev = current2,current2 = (*current2).next);
            (*temp).next = current2;
            if (prev != NULL)
               (*prev).next = temp;
            else
               setvalarray(B,tindex,temp);
         }   
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}





void gaussseidelsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                      int numsteps, double ***S, PBData &pb)
{
// uses values of x as initial guess
   int i, m, n, step, tindex[grid.dim], rindex[grid.dim];
   double temp, diag, ehere;
   SparseElt2 *current;
   
   for (step = 1; step <= numsteps; step++)
   {
      cout << "step = " << step << endl;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
         for (m = 0; m < grid.dim; m++)
            rindex[m] = tindex[m];

         diag = 0.0;
         temp = evalarray(b,tindex);
         if (evalarray(A,tindex) == NULL)
            for (m = 0; m < grid.dim; m++)
            {
               diag += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
               }
               rindex[m] = tindex[m];
            }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               if (sub2ind((*current).cindex,grid.nx,grid.dim) != 
                      sub2ind(tindex,grid.nx,grid.dim))
                  temp -= (*current).val*evalarray(x,(*current).cindex);
               else
                  diag += (*current).val;

         setvalarray(x,tindex,temp/diag);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
}

void gaussseidelsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                      int numsteps, double ***a, double ***S, PBData &pb)
{
// uses values of x as initial guess
   int i, m, n, step, tindex[grid.dim], rindex[grid.dim];
   double temp, diag, ehere;
   SparseElt2 *current;
   
   for (step = 1; step <= numsteps; step++)
   {
      cout << "step = " << step << endl;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
         for (m = 0; m < grid.dim; m++)
            rindex[m] = tindex[m];

         diag = 0.0;
         temp = evalarray(b,tindex);
         if (evalarray(A,tindex) == NULL)
         {
            diag += evalarray(a,rindex);
            for (m = 0; m < grid.dim; m++)
            {
               diag += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
               }
               rindex[m] = tindex[m];
            }
         }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               if (sub2ind((*current).cindex,grid.nx,grid.dim) != 
                      sub2ind(tindex,grid.nx,grid.dim))
                  temp -= (*current).val*evalarray(x,(*current).cindex);
               else
                  diag += (*current).val;

         setvalarray(x,tindex,temp/diag);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
}
//https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
void BICGSTAB(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
              int numsteps, double tol, TempStruct &tmp)
{
   int i, step;
   int tindex[grid.dim];
   double residual, oldresidual, bnorm, temp; 
   double alphanumer, alphadenom, alpha, omeganumer, omegadenom, omega, beta;
   double ***r = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***rstar = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***s = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***p = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***v = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   SparseElt2 *current;

   bnorm = 0.0;
   residual = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(b,tindex);
      for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
         temp -= (*current).val*evalarray(x,(*current).cindex);
      residual += temp*temp; // residual = (b - A 1).^2
      setvalarray(r,tindex,temp); //r = b - A 1
      setvalarray(rstar,tindex,evalarray(r,tindex));// rstar = r = (wiki) r0hat
      setvalarray(p,tindex,evalarray(rstar,tindex));// p1 = r0 
      bnorm += evalarray(b,tindex)*evalarray(b,tindex); // norm of b

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);
   residual = sqrt(residual)/bnorm;
   cout << "   residual = " << residual << endl;
   for (step = 1; step <= numsteps && residual >= tol && 
                  (step == 1 || fabs(oldresidual-residual) > 0.0); step++)
   {
      alphanumer = 0.0;
      alphadenom = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         temp = 0.0;
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
            temp += (*current).val*evalarray(p,(*current).cindex); // temp =  A pi = vi
         alphanumer += evalarray(r,tindex)*evalarray(rstar,tindex); // (r, r0hat)
         alphadenom += temp*evalarray(rstar,tindex); // (vi, r0hat)
         setvalarray(v,tindex,temp);// v = vi

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      alpha = alphanumer/alphadenom;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(s,tindex,evalarray(r,tindex)-alpha*evalarray(v,tindex)); // s = r_(i-1) - alpha v_i

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      omeganumer = 0.0;
      omegadenom = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         temp = 0.0;
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
            temp += (*current).val*evalarray(s,(*current).cindex); // temp = t = A s
         omeganumer += temp*evalarray(s,tindex); // (t,A)
         omegadenom += temp*temp; // (t, t)
         setvalarray(r,tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      omega = omeganumer/omegadenom;

      oldresidual = residual;
      residual = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(x,tindex,evalarray(x,tindex)+alpha*evalarray(p,tindex)+
                                                  omega*evalarray(s,tindex)); // x_i = x_i-1 + alpha p_i + omega_i s
         setvalarray(r,tindex,evalarray(s,tindex)-omega*evalarray(r,tindex)); // r = s - omega_i t
         residual += evalarray(r,tindex)*evalarray(r,tindex); 

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      residual = sqrt(residual)/bnorm;
      if (step%1000 == 0)
         cout << "   step = " << step << " and residual = " << residual 
              << " and change " << oldresidual-residual << " and factor " 
              << oldresidual/residual << endl;

      beta = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         beta += evalarray(r,tindex)*evalarray(rstar,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      beta = beta/alphanumer*alpha/omega;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(p,tindex,evalarray(r,tindex)+
                              beta*(evalarray(p,tindex)-omega*evalarray(v,tindex)));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
   cout << "   step = " << step << " and residual = " << residual << " and change " 
        << oldresidual-residual << " and factor " << oldresidual/residual << endl;

   removefourd(r,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(rstar,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(s,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(p,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(v,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}
// In BICGSTABsmall, at interior point, the coefficients are not added into the sparsematrix
// only difference happens when doing matrix multiplication Ax,
// At an interior point index, A(index) = NULL, need to calculate A(index) x directly without looking for next elemeent
void BICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                   int numsteps, double tol, double ***S, PBData &pb, TempStruct &tmp)
{
   int i, m, n, step;
   int tindex[grid.dim], rindex[grid.dim];
   double residual, oldresidual, bnorm, temp, ehere; 
   double alphanumer, alphadenom, alpha, omeganumer, omegadenom, omega, beta;
   double ***r = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***rstar = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***s = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***p = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***v = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   SparseElt2 *current;
   // ofstream alpha2d("alpha2d.dat", ios::out);

   bnorm = 0.0;
   residual = 0.0;
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

      temp = evalarray(b,tindex);
      if (evalarray(A,tindex) == NULL) // if interior point, A(tindex) x can be calulated directly
         for (m = 0; m < grid.dim; m++)
         {
            temp -= 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            for (n = -1; n <= 1; n += 2)
            {
               rindex[m] = tindex[m]+n;
               if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                  temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            }
            rindex[m] = tindex[m];
         }
      else
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
            temp -= (*current).val*evalarray(x,(*current).cindex);
      residual += temp*temp;
      setvalarray(r,tindex,temp);
      setvalarray(rstar,tindex,evalarray(r,tindex));
      setvalarray(p,tindex,evalarray(rstar,tindex));
      bnorm += evalarray(b,tindex)*evalarray(b,tindex);
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);
   residual = sqrt(residual)/bnorm;
   cout << "   residual = " << residual << endl;
   // alpha2d << residual << endl;
   for (step = 1; step <= numsteps && residual >= tol && 
                  (step == 1 || fabs(oldresidual-residual) > 0.0); step++)
   {
      alphanumer = 0.0;
      alphadenom = 0.0;
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
               temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(p,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(p,rindex);
               }
               rindex[m] = tindex[m];
            }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp += (*current).val*evalarray(p,(*current).cindex);
         alphanumer += evalarray(r,tindex)*evalarray(rstar,tindex);
         alphadenom += temp*evalarray(rstar,tindex);
         setvalarray(v,tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      alpha = alphanumer/alphadenom;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(s,tindex,evalarray(r,tindex)-alpha*evalarray(v,tindex));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      omeganumer = 0.0;
      omegadenom = 0.0;
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
               temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(s,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(s,rindex);
               }
               rindex[m] = tindex[m];
            }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp += (*current).val*evalarray(s,(*current).cindex);
         omeganumer += temp*evalarray(s,tindex);
         omegadenom += temp*temp;
         setvalarray(r,tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      omega = omeganumer/omegadenom;

      oldresidual = residual;
      residual = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(x,tindex,evalarray(x,tindex)+alpha*evalarray(p,tindex)+
                                                  omega*evalarray(s,tindex));
         setvalarray(r,tindex,evalarray(s,tindex)-omega*evalarray(r,tindex));
         residual += evalarray(r,tindex)*evalarray(r,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      residual = sqrt(residual)/bnorm;
      //alpha2d << residual << endl;
      if (step%1000 == 0)
         cout << "   step = " << step << " and residual = " << residual 
              << " and change " << oldresidual-residual << " and factor " 
              << oldresidual/residual << endl;

      beta = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         beta += evalarray(r,tindex)*evalarray(rstar,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      beta = beta/alphanumer*alpha/omega;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(p,tindex,evalarray(r,tindex)+
                              beta*(evalarray(p,tindex)-omega*evalarray(v,tindex)));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
   cout << "   step = " << step << " and residual = " << residual << " and change " 
        << oldresidual-residual << " and factor " << oldresidual/residual << endl;

   removefourd(r,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(rstar,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(s,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(p,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(v,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}

// with a
void BICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                   int numsteps, double tol, double ***a, double ***S, PBData &pb, 
                   TempStruct &tmp)
{
   int i, m, n, step;
   int tindex[grid.dim], rindex[grid.dim];
   double residual, oldresidual, bnorm, temp, ehere; 
   double alphanumer, alphadenom, alpha, omeganumer, omegadenom, omega, beta;
   double ***r = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***rstar = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***s = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***p = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***v = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   SparseElt2 *current;

   bnorm = 0.0;
   residual = 0.0;
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

      temp = evalarray(b,tindex);
      if (evalarray(A,tindex) == NULL)
      {
         temp -= evalarray(a,rindex)*evalarray(x,rindex);
         for (m = 0; m < grid.dim; m++)
         {
            temp -= 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            for (n = -1; n <= 1; n += 2)
            {
               rindex[m] = tindex[m]+n;
               if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                  temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            }
            rindex[m] = tindex[m];
         }
      }
      else
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
            temp -= (*current).val*evalarray(x,(*current).cindex);
      residual += temp*temp;
      setvalarray(r,tindex,temp);
      setvalarray(rstar,tindex,evalarray(r,tindex));
      setvalarray(p,tindex,evalarray(rstar,tindex));
      bnorm += evalarray(b,tindex)*evalarray(b,tindex);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);
   residual = sqrt(residual)/bnorm;
   cout << "   residual = " << residual << endl;
   for (step = 1; step <= numsteps && residual >= tol && 
                  (step == 1 || fabs(oldresidual-residual) > 0.0); step++)
   {
      alphanumer = 0.0;
      alphadenom = 0.0;
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
         {
            temp += evalarray(a,rindex)*evalarray(p,rindex);
            for (m = 0; m < grid.dim; m++)
            {
               temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(p,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(p,rindex);
               }
               rindex[m] = tindex[m];
            }
         }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp += (*current).val*evalarray(p,(*current).cindex);
         alphanumer += evalarray(r,tindex)*evalarray(rstar,tindex);
         alphadenom += temp*evalarray(rstar,tindex);
         setvalarray(v,tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      alpha = alphanumer/alphadenom;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(s,tindex,evalarray(r,tindex)-alpha*evalarray(v,tindex));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      omeganumer = 0.0;
      omegadenom = 0.0;
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
         {
            temp += evalarray(a,rindex)*evalarray(s,rindex);
            for (m = 0; m < grid.dim; m++)
            {
               temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(s,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(s,rindex);
               }
               rindex[m] = tindex[m];
            }
         }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp += (*current).val*evalarray(s,(*current).cindex);
         omeganumer += temp*evalarray(s,tindex);
         omegadenom += temp*temp;
         setvalarray(r,tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      omega = omeganumer/omegadenom;

      oldresidual = residual;
      residual = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(x,tindex,evalarray(x,tindex)+alpha*evalarray(p,tindex)+
                                                  omega*evalarray(s,tindex));
         setvalarray(r,tindex,evalarray(s,tindex)-omega*evalarray(r,tindex));
         residual += evalarray(r,tindex)*evalarray(r,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      residual = sqrt(residual)/bnorm;
      if (step%1000 == 0)
         cout << "   step = " << step << " and residual = " << residual 
              << " and change " << oldresidual-residual << " and factor " 
              << oldresidual/residual << endl;

      beta = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         beta += evalarray(r,tindex)*evalarray(rstar,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      beta = beta/alphanumer*alpha/omega;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(p,tindex,evalarray(r,tindex)+
                              beta*(evalarray(p,tindex)-omega*evalarray(v,tindex)));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
   cout << "   step = " << step << " and residual = " << residual << " and change " 
        << oldresidual-residual << " and factor " << oldresidual/residual << endl;

   removefourd(r,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(rstar,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(s,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(p,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(v,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}

void prerightBICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, 
                           SparseElt2**** &M, GridData &grid, int numsteps, double tol, 
                           double ***S, PBData &pb, TempStruct &tmp)
{
   int i, m, n, step;
   int tindex[grid.dim], rindex[grid.dim];
   double residual, oldresidual, bnorm, temp, ehere; 
   double alphanumer, alphadenom, alpha, omeganumer, omegadenom, omega, beta;
   double ***r = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***rstar = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***s = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***p = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***v = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***y = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   SparseElt2 *current;
//   ofstream alpha2d("alpha2d.dat", ios::out);
   clock_t cstart1, cend1, cstart2, cend2, cstart3, cend3;
   double clocktime1 = 0.0, clocktime2 = 0.0, clocktime3 = 0.0;

   bnorm = 0.0;
   residual = 0.0;
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

      temp = evalarray(b,tindex);
      if (evalarray(A,tindex) == NULL)
         for (m = 0; m < grid.dim; m++)
         {
            temp -= 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            for (n = -1; n <= 1; n += 2)
            {
               rindex[m] = tindex[m]+n;
               if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                  temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            }
            rindex[m] = tindex[m];
         }
      else
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
            temp -= (*current).val*evalarray(x,(*current).cindex);
      residual += temp*temp;
      setvalarray(r,tindex,temp);
      setvalarray(rstar,tindex,evalarray(r,tindex));
      setvalarray(p,tindex,evalarray(rstar,tindex));
      bnorm += evalarray(b,tindex)*evalarray(b,tindex);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);
   residual = sqrt(residual)/bnorm;
   cout << "   residual = " << residual << endl;
//   alpha2d << residual << endl;
   for (step = 1; step <= numsteps && residual >= tol && 
                  (step == 1 || fabs(oldresidual-residual) > 0.0); step++)
   {
      cstart2 = clock ();
      alphanumer = 0.0;
      alphadenom = 0.0;
//      leftmultILUinv(y,M,p,grid);
      cstart1 = clock ();
      leftmultILUinv(s,M,p,grid);
      cend1 = clock ();
      clocktime1 += (double) (cend1-cstart1)/CLOCKS_PER_SEC;
      cstart3 = clock ();
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
               temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(s,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(s,rindex);
               }
               rindex[m] = tindex[m];
            }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp += (*current).val*evalarray(s,(*current).cindex);
         alphanumer += evalarray(r,tindex)*evalarray(rstar,tindex);
         alphadenom += temp*evalarray(rstar,tindex);
         setvalarray(v,tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      cend3 = clock ();
      clocktime3 += (double) (cend3-cstart3)/CLOCKS_PER_SEC;
      alpha = alphanumer/alphadenom;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(x,tindex,evalarray(x,tindex)+alpha*evalarray(s,tindex));
         setvalarray(s,tindex,evalarray(r,tindex)-alpha*evalarray(v,tindex));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      omeganumer = 0.0;
      omegadenom = 0.0;
//      leftmultILUinv(z,M,s,grid);
      cstart1 = clock ();
      leftmultILUinv(y,M,s,grid);
      cend1 = clock ();
      clocktime1 += (double) (cend1-cstart1)/CLOCKS_PER_SEC;
      cstart3 = clock ();
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
               temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(y,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(y,rindex);
               }
               rindex[m] = tindex[m];
            }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp += (*current).val*evalarray(y,(*current).cindex);
         omeganumer += temp*evalarray(s,tindex);
         omegadenom += temp*temp;
         setvalarray(r,tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      cend3 = clock ();
      clocktime3 += (double) (cend3-cstart3)/CLOCKS_PER_SEC;
      omega = omeganumer/omegadenom;

      oldresidual = residual;
      residual = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(x,tindex,evalarray(x,tindex)+omega*evalarray(y,tindex));
//         setvalarray(x,tindex,evalarray(x,tindex)+alpha*evalarray(y,tindex)+
//                                                  omega*evalarray(z,tindex));
         setvalarray(r,tindex,evalarray(s,tindex)-omega*evalarray(r,tindex));
         residual += evalarray(r,tindex)*evalarray(r,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      residual = sqrt(residual)/bnorm;
//      alpha2d << residual << endl;
      if (step%10 == 0)
         cout << "   step = " << step << " and residual = " << residual 
              << " and change " << oldresidual-residual << " and factor " 
              << oldresidual/residual << endl;

      beta = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         beta += evalarray(r,tindex)*evalarray(rstar,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      beta = beta/alphanumer*alpha/omega;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(p,tindex,evalarray(r,tindex)+
                              beta*(evalarray(p,tindex)-omega*evalarray(v,tindex)));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      cend2 = clock ();
      clocktime2 += (double) (cend2-cstart2)/CLOCKS_PER_SEC;
//      cout << step << " " << residual << endl;
//      getchar();
   }
   cout << "   step = " << step-1 << " and residual = " << residual << " and change " 
        << oldresidual-residual << " and factor " << oldresidual/residual << endl;
   cout << "ILU inv time per = " << clocktime1/(step-1)/2.0 << " out of "
        << clocktime2/(step-1) << " compared to " << clocktime3/(step-1)/2.0 << endl;

   removefourd(r,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(rstar,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(s,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(p,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(v,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(y,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}

void prerightBICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, 
                           SparseElt2**** &M, GridData &grid, int numsteps, double tol, 
                           double ***a, double ***S, PBData &pb, TempStruct &tmp)
{
   int i, m, n, step;
   int tindex[grid.dim], rindex[grid.dim];
   double residual, oldresidual, bnorm, temp, ehere; 
   double alphanumer, alphadenom, alpha, omeganumer, omegadenom, omega, beta;
   double ***r = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***rstar = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***s = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***p = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***v = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***y = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   SparseElt2 *current;
//   ofstream alpha2d("alpha2d.dat", ios::out);
   clock_t cstart1, cend1, cstart2, cend2, cstart3, cend3;
   double clocktime1 = 0.0, clocktime2 = 0.0, clocktime3 = 0.0;

   bnorm = 0.0;
   residual = 0.0;
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

      temp = evalarray(b,tindex);
      if (evalarray(A,tindex) == NULL)
      {
         temp -= evalarray(a,rindex)*evalarray(x,rindex);
         for (m = 0; m < grid.dim; m++)
         {
            temp -= 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            for (n = -1; n <= 1; n += 2)
            {
               rindex[m] = tindex[m]+n;
               if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                  temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            }
            rindex[m] = tindex[m];
         }
      }
      else
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
            temp -= (*current).val*evalarray(x,(*current).cindex);
      residual += temp*temp;
      setvalarray(r,tindex,temp);
      setvalarray(rstar,tindex,evalarray(r,tindex));
      setvalarray(p,tindex,evalarray(rstar,tindex));
      bnorm += evalarray(b,tindex)*evalarray(b,tindex);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);
   residual = sqrt(residual)/bnorm;
   cout << "   residual = " << residual << endl;
//   alpha2d << residual << endl;
   for (step = 1; step <= numsteps && residual >= tol && 
                  (step == 1 || fabs(oldresidual-residual) > 0.0); step++)
   {
      cstart2 = clock ();
      alphanumer = 0.0;
      alphadenom = 0.0;
//      leftmultILUinv(y,M,p,grid);
      cstart1 = clock ();
      leftmultILUinv(s,M,p,grid);
      cend1 = clock ();
      clocktime1 += (double) (cend1-cstart1)/CLOCKS_PER_SEC;
      cstart3 = clock ();
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
         {
            temp += evalarray(a,rindex)*evalarray(s,rindex);
            for (m = 0; m < grid.dim; m++)
            {
               temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(s,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(s,rindex);
               }
               rindex[m] = tindex[m];
            }
         }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp += (*current).val*evalarray(s,(*current).cindex);
         alphanumer += evalarray(r,tindex)*evalarray(rstar,tindex);
         alphadenom += temp*evalarray(rstar,tindex);
         setvalarray(v,tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      cend3 = clock ();
      clocktime3 += (double) (cend3-cstart3)/CLOCKS_PER_SEC;
      alpha = alphanumer/alphadenom;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(x,tindex,evalarray(x,tindex)+alpha*evalarray(s,tindex));
         setvalarray(s,tindex,evalarray(r,tindex)-alpha*evalarray(v,tindex));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      omeganumer = 0.0;
      omegadenom = 0.0;
//      leftmultILUinv(z,M,s,grid);
      cstart1 = clock ();
      leftmultILUinv(y,M,s,grid);
      cend1 = clock ();
      clocktime1 += (double) (cend1-cstart1)/CLOCKS_PER_SEC;
      cstart3 = clock ();
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
         {
            temp += evalarray(a,rindex)*evalarray(y,rindex);
            for (m = 0; m < grid.dim; m++)
            {
               temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(y,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(y,rindex);
               }
               rindex[m] = tindex[m];
            }
         }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp += (*current).val*evalarray(y,(*current).cindex);
         omeganumer += temp*evalarray(s,tindex);
         omegadenom += temp*temp;
         setvalarray(r,tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      cend3 = clock ();
      clocktime3 += (double) (cend3-cstart3)/CLOCKS_PER_SEC;
      omega = omeganumer/omegadenom;

      oldresidual = residual;
      residual = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(x,tindex,evalarray(x,tindex)+omega*evalarray(y,tindex));
//         setvalarray(x,tindex,evalarray(x,tindex)+alpha*evalarray(y,tindex)+
//                                                  omega*evalarray(z,tindex));
         setvalarray(r,tindex,evalarray(s,tindex)-omega*evalarray(r,tindex));
         residual += evalarray(r,tindex)*evalarray(r,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      residual = sqrt(residual)/bnorm;
//      alpha2d << residual << endl;
      if (step%10 == 0)
         cout << "   step = " << step << " and residual = " << residual 
              << " and change " << oldresidual-residual << " and factor " 
              << oldresidual/residual << endl;

      beta = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         beta += evalarray(r,tindex)*evalarray(rstar,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      beta = beta/alphanumer*alpha/omega;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(p,tindex,evalarray(r,tindex)+
                              beta*(evalarray(p,tindex)-omega*evalarray(v,tindex)));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      cend2 = clock ();
      clocktime2 += (double) (cend2-cstart2)/CLOCKS_PER_SEC;
//      cout << step << " " << residual << endl;
//      getchar();
   }
   cout << "   step = " << step-1 << " and residual = " << residual << " and change " 
        << oldresidual-residual << " and factor " << oldresidual/residual << endl;
   cout << "ILU inv time per = " << clocktime1/(step-1)/2.0 << " out of "
        << clocktime2/(step-1) << " compared to " << clocktime3/(step-1)/2.0 << endl;

   removefourd(r,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(rstar,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(s,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(p,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(v,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(y,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}

void preBICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, SparseElt2**** M, 
                      GridData &grid, int numsteps, double tol, double ***S, 
                      PBData &pb, TempStruct &tmp)
{
   int i, m, n, step;
   int tindex[grid.dim], rindex[grid.dim];
   double residual, oldresidual, bnorm, temp, ehere; 
   double alphanumer, alphadenom, alpha, omeganumer, omegadenom, omega, beta;
   double ***r = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***rstar = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***s = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***p = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid),
          ***v = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   SparseElt2 *current;
   // ofstream alpha2d("alpha2d.dat", ios::out);

   leftmultILUinv(p,M,b,grid);
   leftmultmxsmall(rstar,A,x,grid,S,pb);
   leftmultILUinv(r,M,rstar,grid);
   bnorm = 0.0;
   residual = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(r,tindex,evalarray(p,tindex)-evalarray(r,tindex));
      bnorm += evalarray(p,tindex)*evalarray(p,tindex);
      residual += evalarray(r,tindex)*evalarray(r,tindex);
      setvalarray(rstar,tindex,evalarray(r,tindex));
      setvalarray(p,tindex,evalarray(rstar,tindex));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);
   residual = sqrt(residual)/bnorm;
   cout << "   residual = " << residual << endl;
   //alpha2d << residual << endl;
   for (step = 1; step <= numsteps && residual >= tol && 
                  (step == 1 || fabs(oldresidual-residual) > 0.0); step++)
   {
      leftmultmxsmall(v,A,p,grid,S,pb);
      leftmultILUinv(v,M,v,grid);
      alphanumer = 0.0;
      alphadenom = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         alphanumer += evalarray(r,tindex)*evalarray(rstar,tindex);
         alphadenom += evalarray(v,tindex)*evalarray(rstar,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      alpha = alphanumer/alphadenom;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(s,tindex,evalarray(r,tindex)-alpha*evalarray(v,tindex));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      omeganumer = 0.0;
      omegadenom = 0.0;
      leftmultmxsmall(r,A,s,grid,S,pb);
      leftmultILUinv(r,M,r,grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         omeganumer += evalarray(r,tindex)*evalarray(s,tindex);
         omegadenom += evalarray(r,tindex)*evalarray(r,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      omega = omeganumer/omegadenom;

      oldresidual = residual;
      residual = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(x,tindex,evalarray(x,tindex)+alpha*evalarray(p,tindex)+
                                                  omega*evalarray(s,tindex));
         setvalarray(r,tindex,evalarray(s,tindex)-omega*evalarray(r,tindex));
         residual += evalarray(r,tindex)*evalarray(r,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      residual = sqrt(residual)/bnorm;
      //alpha2d << residual << endl;
      if (step%1000 == 0)
         cout << "   step = " << step << " and residual = " << residual 
              << " and change " << oldresidual-residual << " and factor " 
              << oldresidual/residual << endl;

      beta = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         beta += evalarray(r,tindex)*evalarray(rstar,tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      beta = beta/alphanumer*alpha/omega;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(p,tindex,evalarray(r,tindex)+
                              beta*(evalarray(p,tindex)-omega*evalarray(v,tindex)));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
   cout << "   step = " << step << " and residual = " << residual << " and change " 
        << oldresidual-residual << " and factor " << oldresidual/residual << endl;

   removefourd(r,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(rstar,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(s,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(p,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(v,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}




void ILU(SparseElt2**** &M, SparseElt2**** &A, GridData &grid)
// should write an ILUsmall also
{
   int i, r, tindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   int N = 1;
   for (r = 0; r < grid.dim; r++)
      N *= grid.nx[r]+1;

// copy A to M and order M
   clearsparse(M,grid);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(M,tindex,NULL);
      for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
      {
         temp = new SparseElt2;
         (*temp).cindex = new int[grid.dim];
         for (r = 0; r < grid.dim; r++) 
            (*temp).cindex[r] = (*current).cindex[r];
         (*temp).val = (*current).val;
         for (prev = NULL,current2 = evalarray(M,tindex); 
              current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                  sub2ind((*temp).cindex,grid.nx,grid.dim); 
              prev = current2,current2 = (*current2).next);
         (*temp).next = current2;
         if (prev != NULL)
            (*prev).next = temp;
         else
            setvalarray(M,tindex,temp);
      }   
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   
// work on M
   for (i = 0; i < N; i++)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
// row is tindex, column is (*current).cindex
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
      {
// (*current2).cindex is different columns of row tindex
         for (current2 = evalarray(M,tindex); current2 != current && 
              sub2ind((*current2).cindex,grid.nx,grid.dim) < 
              min(sub2ind(tindex,grid.nx,grid.dim),
                  sub2ind((*current).cindex,grid.nx,grid.dim)); 
              current2 = (*current2).next)
         {
// (*temp).cindex is different columns of row (*current2).cindex
            for (temp = evalarray(M,(*current2).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val -= ((*current2).val)*((*temp).val);
         }
// if working on L
         if (sub2ind((*current).cindex,grid.nx,grid.dim) < 
             sub2ind(tindex,grid.nx,grid.dim))
         {
            for (temp = evalarray(M,(*current).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val /= (*temp).val;
         }
      }
   }
}




void ILUsmall(SparseElt2**** &M, SparseElt2**** &A, GridData &grid, double ***S, 
              PBData &pb)
{
   int i, r, m, n, tindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   double ehere;
   int N = 1;
   for (r = 0; r < grid.dim; r++)
      N *= grid.nx[r]+1;

   copyfromsmall(M,A,grid,S,pb);
/*
// copy A to M and order M
   clearsparse(M,grid);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(M,tindex,NULL);
      if (evalarray(A,tindex) == NULL)
      {
// create rows in M when they don't exist in A
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
         temp = new SparseElt2;
         (*temp).cindex = new int[grid.dim];
         for (r = 0; r < grid.dim; r++) 
            (*temp).cindex[r] = tindex[r];
         (*temp).val = 0.0;
         for (m = 0; m < grid.dim; m++)
            (*temp).val += 2.0*ehere/(grid.dx[m]*grid.dx[m]);
         (*temp).next = NULL;
         setvalarray(M,tindex,temp);

         for (m = 0; m < grid.dim; m++)
            for (n = -1; n <= 1; n += 2)
            {
               temp = new SparseElt2;
               (*temp).cindex = new int[grid.dim];
               for (r = 0; r < grid.dim; r++) 
                  (*temp).cindex[r] = tindex[r];
               (*temp).cindex[m] = tindex[m]+n;
               if ((*temp).cindex[m] >= 0 && (*temp).cindex[m] <= grid.nx[m])
                  (*temp).val = -ehere/(grid.dx[m]*grid.dx[m]);
               for (prev = NULL,current2 = evalarray(M,tindex); 
                    current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                        sub2ind((*temp).cindex,grid.nx,grid.dim); 
                    prev = current2,current2 = (*current2).next);
               (*temp).next = current2;
               if (prev != NULL)
                  (*prev).next = temp;
               else
                  setvalarray(M,tindex,temp);
            }
      }
      else
         for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
         {
            temp = new SparseElt2;
            (*temp).cindex = new int[grid.dim];
            for (r = 0; r < grid.dim; r++) 
               (*temp).cindex[r] = (*current).cindex[r];
            (*temp).val = (*current).val;
            for (prev = NULL,current2 = evalarray(M,tindex); 
                 current2 != NULL && sub2ind((*current2).cindex,grid.nx,grid.dim) < 
                                     sub2ind((*temp).cindex,grid.nx,grid.dim); 
                 prev = current2,current2 = (*current2).next);
            (*temp).next = current2;
            if (prev != NULL)
               (*prev).next = temp;
            else
               setvalarray(M,tindex,temp);
         }   
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
*/

// work on M
   for (i = 0; i < N; i++)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
// row is tindex, column is (*current).cindex
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
      {
// (*current2).cindex is different columns of row tindex
         for (current2 = evalarray(M,tindex); current2 != current && 
              sub2ind((*current2).cindex,grid.nx,grid.dim) < 
              min(sub2ind(tindex,grid.nx,grid.dim),
                  sub2ind((*current).cindex,grid.nx,grid.dim)); 
              current2 = (*current2).next)
         {
// (*temp).cindex is different columns of row (*current2).cindex
            for (temp = evalarray(M,(*current2).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val -= ((*current2).val)*((*temp).val);
         }
// if working on L
         if (sub2ind((*current).cindex,grid.nx,grid.dim) < 
             sub2ind(tindex,grid.nx,grid.dim))
         {
            for (temp = evalarray(M,(*current).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val /= (*temp).val;
         }
      }
   }
}

void ILUsmall(SparseElt2**** &M, SparseElt2**** &A, GridData &grid, double ***a, 
              double ***S, PBData &pb)
{
   int i, r, m, n, tindex[grid.dim];
   SparseElt2 *current, *current2, *prev, *temp;
   double ehere;
   int N = 1;
   for (r = 0; r < grid.dim; r++)
      N *= grid.nx[r]+1;

   copyfromsmall(M,A,grid,a,S,pb);

// work on M
   for (i = 0; i < N; i++)
   {
      ind2sub(tindex,i,grid.nx,grid.dim);
// row is tindex, column is (*current).cindex
      for (current = evalarray(M,tindex); current != NULL; current = (*current).next)
      {
// (*current2).cindex is different columns of row tindex
         for (current2 = evalarray(M,tindex); current2 != current && 
              sub2ind((*current2).cindex,grid.nx,grid.dim) < 
              min(sub2ind(tindex,grid.nx,grid.dim),
                  sub2ind((*current).cindex,grid.nx,grid.dim)); 
              current2 = (*current2).next)
         {
// (*temp).cindex is different columns of row (*current2).cindex
            for (temp = evalarray(M,(*current2).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val -= ((*current2).val)*((*temp).val);
         }
// if working on L
         if (sub2ind((*current).cindex,grid.nx,grid.dim) < 
             sub2ind(tindex,grid.nx,grid.dim))
         {
            for (temp = evalarray(M,(*current).cindex); temp != NULL && 
                 sub2ind((*temp).cindex,grid.nx,grid.dim) < 
                 sub2ind((*current).cindex,grid.nx,grid.dim); temp = (*temp).next);
            if (temp != NULL && sub2ind((*temp).cindex,grid.nx,grid.dim) == 
                                sub2ind((*current).cindex,grid.nx,grid.dim))
               (*current).val /= (*temp).val;
         }
      }
   }
}

