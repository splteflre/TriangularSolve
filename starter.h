#include <iostream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm> 
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <omp.h>
#include <math.h>
#include <cfloat>


void sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double **x);
void parallel_sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x, int nlev, int *ilev, int *jlev);
void spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax, const double *x, double *y);
void create_csc(char *matrix, char *b, int **Lp, int **Li, double **Lx, int &n, double **x, std::vector<int> &non_zero);
void create_level_set(int n, int **Lp, int **Li, int **jlev, int **ilev, int &nlev, std::vector<int> non_zero);
bool AreSame(double a, double b, double epsilon);
