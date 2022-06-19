#include "advance.h"
#include "tryvn.h"

extern int globrk;
extern int globexactmotion;
extern int globexactvn;
extern int globperturb;

using namespace std;

//u is surface, update u with temp
void advance(double ***u, PBData &pb, MarchStruct &march, TempStruct &tmp, 
             GridData &grid)
{
   int i, tindex[grid.dim];
   double ***temp = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);

   if (!globrk)
   {
      if (!globexactmotion && !globexactvn)//not exact motion and vn
         getRHS(temp,u,pb,march,tmp,grid);
      else if (globexactmotion && !globexactvn)
         getRHShyp(temp,u,pb,march,tmp,grid);
      else if (globexactvn)
         getRHSexactvn(temp,u,pb,march,tmp,grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(u,tindex,evalarray(u,tindex)+grid.dt*evalarray(temp,tindex));
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      grid.t += grid.dt;
      grid.radius = getexactradius(grid.t,grid.radius0,grid.radius,1.0e-14,100,grid);
   }
   else
   {
      double tempr, tempt, tempdt;

      if (!globexactmotion && !globexactvn)
         getRHS(temp,u,pb,march,tmp,grid);
      else if (globexactmotion && !globexactvn)
         getRHShyp(temp,u,pb,march,tmp,grid);
      else if (globexactvn)
         getRHSexactvn(temp,u,pb,march,tmp,grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(temp,tindex,evalarray(u,tindex)+
                                 grid.dt*evalarray(temp,tindex));
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      if (globperturb == 0)
         perturb(temp,grid.tol,pb,grid);
      else if (globperturb > 0)
         perturbstatus(temp,grid.tol,globperturb,pb,grid);
      cout << "max vn error = " << grid.maxvn << endl;
      //delta2d << grid.t << " " << grid.maxvn << endl;

      double ***temp2 = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);

      tempt = grid.t;
      tempr = grid.radius;
      tempdt = grid.dt;

      grid.t += grid.dt;
      grid.radius = getexactradius(grid.t,grid.radius0,grid.radius,1.0e-14,100,grid);
      if (!globexactmotion && !globexactvn)
         getRHS(temp2,temp,pb,march,tmp,grid);
      else if (globexactmotion && !globexactvn)
         getRHShyp(temp2,temp,pb,march,tmp,grid);
      else if (globexactvn)
         getRHSexactvn(temp2,temp,pb,march,tmp,grid);

//      grid.t = tempt;
//      grid.radius = tempr;
      grid.dt = tempdt;
      cout << "   actually using in 2nd step dt = " << grid.dt << endl;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(u,tindex,0.5*evalarray(u,tindex)+0.5*evalarray(temp,tindex)+
                              0.5*grid.dt*evalarray(temp2,tindex));
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      removefourd(temp2,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   }

   removefourd(temp,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}

void advance(double ***u, double ***a, PBData &pb, MarchStruct &march, TempStruct &tmp, 
             GridData &grid)
{
   int i, tindex[grid.dim];
   double ***temp = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);

   if (!globrk)
   {
      if (!globexactmotion && !globexactvn)
         getRHS(temp,u,a,pb,march,tmp,grid);
      else if (globexactmotion && !globexactvn)
         getRHShyp(temp,u,pb,march,tmp,grid);
      else if (globexactvn)
         getRHSexactvn(temp,u,pb,march,tmp,grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(u,tindex,evalarray(u,tindex)+grid.dt*evalarray(temp,tindex));
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      grid.t += grid.dt;
      grid.radius = getexactradius(grid.t,grid.radius0,grid.radius,1.0e-14,100,grid);
   }
   else
   {
      double tempr, tempt, tempdt;

      if (!globexactmotion && !globexactvn)
         getRHS(temp,u,pb,march,tmp,grid);
      else if (globexactmotion && !globexactvn)
         getRHShyp(temp,u,pb,march,tmp,grid);
      else if (globexactvn)
         getRHSexactvn(temp,u,pb,march,tmp,grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(temp,tindex,evalarray(u,tindex)+
                                 grid.dt*evalarray(temp,tindex));
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      if (globperturb == 0)
         perturb(temp,grid.tol,pb,grid);
      else if (globperturb > 0)
         perturbstatus(temp,grid.tol,globperturb,pb,grid);
      cout << "max vn error = " << grid.maxvn << endl;
      //delta2d << grid.t << " " << grid.maxvn << endl;

      double ***temp2 = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);

      tempt = grid.t;
      tempr = grid.radius;
      tempdt = grid.dt;

      grid.t += grid.dt;
      grid.radius = getexactradius(grid.t,grid.radius0,grid.radius,1.0e-14,100,grid);
      if (!globexactmotion && !globexactvn)
         getRHS(temp2,temp,pb,march,tmp,grid);
      else if (globexactmotion && !globexactvn)
         getRHShyp(temp2,temp,pb,march,tmp,grid);
      else if (globexactvn)
         getRHSexactvn(temp2,temp,pb,march,tmp,grid);

//      grid.t = tempt;
//      grid.radius = tempr;
      grid.dt = tempdt;
      cout << "   actually using in 2nd step dt = " << grid.dt << endl;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(u,tindex,0.5*evalarray(u,tindex)+0.5*evalarray(temp,tindex)+
                              0.5*grid.dt*evalarray(temp2,tindex));
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      removefourd(temp2,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   }

   removefourd(temp,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}


//u is extended velocity, one euler step of du/dt = laplacian(u), u(t+1) = u(t) + dt(uxx+uyy+uzz)
void advanceheat(double ***u, TempStruct &tmp, GridData &grid)
{
   int i, r, s, tindex[grid.dim], rindex[grid.dim];
   double deriv;
   double dtheat = grid.mindx*grid.mindx/(2.0*grid.dim); //dt = dx^2/(2*dim)/2
   double ***temp = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(temp,tindex,0.0);
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         deriv = 0.0;
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = tindex[r]+s;
            if (rindex[r] >= 0 && rindex[r] <= grid.nx[r])
               deriv += evalarray(u,rindex);
            else
               deriv += evalarray(u,tindex); //treat boundary as current, von Neumann BC, gradient = 0
         }
         rindex[r] = tindex[r];
         deriv += -2.0*evalarray(u,rindex);
         deriv /= grid.dx[r]*grid.dx[r]; //deriv is second derivative by central differencing
         setvalarray(temp,tindex,evalarray(temp,tindex)+deriv);
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(u,tindex,evalarray(u,tindex)+dtheat*evalarray(temp,tindex));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cout<<"heatsmooth with dt = " <<dtheat<<endl;
   removefourd(temp,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}
//u is extended velocity, solve du/dt = laplacian(u) unitl finaltime, u(t+1) = u(t) + dt(uxx+uyy+uzz)
void advanceheat(double ***u, double finaltime, TempStruct &tmp, GridData &grid)
{
   int i, r, s, tindex[grid.dim], rindex[grid.dim];
   double deriv, dtheat, thetime;
   double ***temp = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);

   thetime = 0.0;
   while (thetime < finaltime)
   {
      dtheat = grid.mindx*grid.mindx/(2.0*grid.dim);
      if (thetime+dtheat > finaltime)
         dtheat = finaltime-thetime;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(temp,tindex,0.0);
         for (r = 0; r < grid.dim; r++)
            rindex[r] = tindex[r];
         for (r = 0; r < grid.dim; r++)
         {
            deriv = 0.0;
            for (s = -1; s <= 1; s += 2)
            {
               rindex[r] = tindex[r]+s;
               if (rindex[r] >= 0 && rindex[r] <= grid.nx[r])
                  deriv += evalarray(u,rindex);
               else
                  deriv += evalarray(u,tindex);
            }
            rindex[r] = tindex[r];
            deriv += -2.0*evalarray(u,rindex);
            deriv /= grid.dx[r]*grid.dx[r];
            setvalarray(temp,tindex,evalarray(temp,tindex)+deriv);
         }

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(u,tindex,evalarray(u,tindex)+dtheat*evalarray(temp,tindex));
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      cout<<"heatsmooth with dt = " <<dtheat<<endl;

      thetime += dtheat;
   }
   removefourd(temp,tmp.fourd,tmp.fdstatus,tmp.Nfd);
}


