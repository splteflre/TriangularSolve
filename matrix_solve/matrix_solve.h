#include <omp.h>
#include <stdio.h>
#include <algorithm> 
#include <vector>
#include <math.h>

void sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double **x);
void parallel_sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x, int nlev, int *ilev, int *jlev);
void spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax, const double *x, double *y);
