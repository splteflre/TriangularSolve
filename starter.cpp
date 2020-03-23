#include "starter.h"

/*
 * Default double comparison is not accurate enough using the code comparison found here instead:
 * https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
 */
bool AreSame(double a, double b, double epsilon)
{
    // Testing different floating point comparison methods
    //return a == b;
    //double relth = FLT_MIN;
    //auto diff = abs(a-b);
    //auto norm = min((abs(a) + abs(b)), numeric_limits<double>::max());
    //return diff < max(relth, epsilon * norm);

    return fabs(a - b) < epsilon;
}

int main(int argc, char *argv[]){
    
    printf("argc: %d\n", argc);
    if(argc != 3){
        printf("wrong # of args\n");
        return 0;
    }

    //TODO: read a matrix from Suite sparse matrix repo
    double *Lx; 
    int *Lp, *Li;                           // Lx, Lp, and Li arrays as described in the handout
    double *x;                              // Right hand side vector, b
    int n;                                // Dimension of the matrix read from input file
    int *jlev;                              // Array that lists the unknowns in ascending order of their levels
    int *ilev;                              // Array that contains the pointers to the levels in jlev
    int nlev;                               // nlev is the number of levels
    double epsilon = 128 * FLT_EPSILON;     // Epsilon used for comparing flaots
    double *tmp;                            // Array used for sanity checking

    std::vector<int> non_zero;
    create_csc(argv[1], argv[2], &Lp, &Li, &Lx, n, &x, non_zero);
    create_level_set(n, &Lp, &Li, &jlev, &ilev, nlev, non_zero);

    //Copying right handside vector for sanity check
    double *y = new double[n];              // vector used for sanity check
    for(int i = 0; i < n; i++){
       y[i] = x[i]; 
    }

    //TODO: Does a triangular solve
    //sptrsv_csc(n, Lp, Li, Lx, &x);
    parallel_sptrsv_csc(n, Lp, Li, Lx, x, nlev, ilev, jlev);
    

    //TODO: sanity check using spmv
    tmp = new double[n];
    std::fill(tmp, tmp+n, 0.);
    spmv_csc(n, Lp, Li, Lx, x, tmp);

    
    for(int i = 0; i < n; i++){
        if (!AreSame(y[i], tmp[i], epsilon)){
            printf("failed for unknown %d\n", i+1);
            printf("x: %f, y: %f\n", y[i], tmp[i]);
            return 0;
        }
    }
    printf("Passed!\n");

    delete[] Lx;
    delete[] Lp;
    delete[] Li;
    delete[] x;
    delete[] tmp;


    return 1;
}
