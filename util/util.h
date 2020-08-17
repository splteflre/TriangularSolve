#include <cfloat>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <string.h>

namespace util 
{
    bool AreSame(double a, double b, double epsilon);
    void create_level_set(int n, int **Lp, int **Li, int **jlev, int **ilev, int &nlev, std::vector<int> non_zero, char* format);
    void init_arrays(char *matrix, char *b, int **Lp, int **Li, double **Lx, int &n, double **x, std::vector<int> &non_zero, char *format, double **d);
};
