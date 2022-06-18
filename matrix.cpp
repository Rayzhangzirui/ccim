#include "matrix.h"
#include <iostream>
using namespace std;

double **matrix(int row, int col)
{
   int i;
   double **m;
   int err = 0;

   m = new double*[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new double[col+1]();
      if (!m[i] && !err)
         err = 1;
   }
   if (err)
      cout << "matrix stuff didn`t work" << endl;
   return m;
}

double ***matrix(int row, int col, int fr)
{
   int i, j;
   double ***m;
   int err = 0;

   m = new double**[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new double*[col+1];
      if (!m[i] && !err)
         err = 1;
      for (j = 0; j <= col; j++)
      {
         m[i][j] = new double[fr+1]();
         if (!m[i][j] && !err)
            err = 1;
      }
   }
   if (err)
      cout << "matrix stuff didn`t work" << endl;
   return m;
}

double ****matrix(int row, int col, int fr, int ss)
{
   int i, j, k;
   double ****m;
   int err = 0;

   m = new double***[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new double**[col+1];
      if (!m[i] && !err)
         err = 1;
      for (j = 0; j <= col; j++)
      {
         m[i][j] = new double*[fr+1];
         if (!m[i][j] && !err)
            err = 1;
         for (k = 0; k <= fr; k++)
         {
            m[i][j][k] = new double[ss+1]();
            if (!m[i][j][k] && !err)
               err = 1;
         }
      }
   }
   if (err)
      cout << "matrix stuff didn`t work" << endl;
   return m;
}

char **cmatrix(int row, int col)
{
   int i;
   char **m;
   int err = 0;

   m = new char*[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new char[col+1];
      if (!m[i] && !err)
         err = 1;
   }
   if (err)
      cout << "matrix stuff didn`t work" << endl;
   return m;
}

char ***cmatrix(int row, int col, int fr)
{
   int i, j;
   char ***m;
   int err = 0;

   m = new char**[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new char*[col+1];
      if (!m[i] && !err)
         err = 1;
      for (j = 0; j <= col; j++)
      {
         m[i][j] = new char[fr+1];
         if (!m[i][j] && !err)
            err = 1;
      }
   }
   if (err)
      cout << "matrix stuff didn`t work" << endl;
   return m;
}

int **imatrix(int row, int col)
{
   int i;
   int **m;
   int err = 0;

   m = new int*[row+1];
   if (!m)
      err = 1;
   for (i = 0; i <= row; i++)
   {
      m[i] = new int[col+1];
      if (!m[i] && !err)
         err = 1;
   }
   if (err)
      cout << "matrix stuff didn`t work" << endl;
   return m;
}


void free_matrix(double **m, int row, int col)
{
   int i;

   for (i = 0; i <= row; i++)
      delete [] m[i];
   delete [] m;
}

void free_matrix(int **m, int row, int col)
{
   int i;

   for (i = 0; i <= row; i++)
      delete [] m[i];
   delete [] m;
}

void free_matrix(char **m, int row, int col)
{
   int i;

   for (i = 0; i <= row; i++)
      delete [] m[i];
   delete [] m;
}

void free_matrix(char ***m, int row, int col, int frb)
{
   int i, j;

   for (i = 0; i <= row; i++)
   {
      for (j = 0; j <= col; j++)
         delete [] m[i][j];
      delete [] m[i];
   }
   delete [] m;
}

void free_matrix(double ***m, int row, int col, int frb)
{
   int i, j;

   for (i = 0; i <= row; i++)
   {
      for (j = 0; j <= col; j++)
         delete [] m[i][j];
      delete [] m[i];
   }
   delete [] m;
}

void free_matrix(double ****m, int row, int col, int fr, int sb)
{
   int i, j, k, l;

   for (i = 0; i <= row; i++)
   {
      for (j = 0; j <= col; j++)
      {
         for (k = 0; k <= fr; k++)
            delete [] m[i][j][k];
         delete [] m[i][j];
      }
      delete [] m[i];
   }
   delete [] m;
}