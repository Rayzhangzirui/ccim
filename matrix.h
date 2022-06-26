#ifndef MATRIX_H
#define MATRIX_H

double **matrix(int row, int col);
double ***matrix(int row, int col, int fr);
double ****matrix(int row, int col, int fr, int ss);
char **cmatrix(int row, int col);
char ***cmatrix(int row, int col, int fr);
int **imatrix(int row, int col);

void free_matrix(double **m, int row, int col);
void free_matrix(int **m, int row, int col);
void free_matrix(char **m, int row, int col);
void free_matrix(char ***m, int row, int col, int frb);
void free_matrix(double ***m, int row, int col, int frb);
void free_matrix(double ****m, int row, int col, int fr, int sb);
inline double evalarray(double ****A, int *index);
inline double evalarray(double ***A, int *index);
inline char evalarray(char ***A, int *index);
inline void setvalarray(char ***A, int *index, char value);
inline double evalarray(double **A, int *index);
inline char evalarray(char **A, int *index);
inline int evalarray(int ***A, int *index);
inline void setvalarray(int ***A, int *index, int value);
inline void setvalarray(char **A, int *index, char value);
inline void setvalarray(double ****A, int *index, double value);
inline void setvalarray(double ***A, int *index, double value);
inline void setvalarray(double **A, int *index, double value);


inline double evalarray(double ****A, int *index)
{
   return A[index[0]][index[1]][index[2]][index[3]];
}

inline double evalarray(double ***A, int *index)
{
   return A[index[0]][index[1]][index[2]];
}

inline char evalarray(char ***A, int *index)
{
   return A[index[0]][index[1]][index[2]];
}

inline void setvalarray(char ***A, int *index, char value)
{
   A[index[0]][index[1]][index[2]] = value;
}

inline double evalarray(double **A, int *index)
{
   return A[index[0]][index[1]];
}

inline char evalarray(char **A, int *index)
{
   return A[index[0]][index[1]];
}

inline int evalarray(int ***A, int *index)
{
   return A[index[0]][index[1]][index[2]];
}

inline void setvalarray(int ***A, int *index, int value)
{
   A[index[0]][index[1]][index[2]] = value;
}



inline void setvalarray(char **A, int *index, char value)
{
   A[index[0]][index[1]] = value;
}

inline void setvalarray(double ****A, int *index, double value)
{
   A[index[0]][index[1]][index[2]][index[3]] = value;
}

inline void setvalarray(double ***A, int *index, double value)
{
   A[index[0]][index[1]][index[2]] = value;
}

inline void setvalarray(double **A, int *index, double value)
{
   A[index[0]][index[1]] = value;
}

inline int equal(int a[], int b[]){
   return a[0]==b[0] && a[1]==b[1] && a[2]==b[2];
}

#endif