#include "util/util.h"
#include "matrix_solve/matrix_solve.h"


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
    int n;                                  // Dimension of the matrix read from input file
    int *jlev;                              // Array that lists the unknowns in ascending order of their levels
    int *ilev;                              // Array that contains the pointers to the levels in jlev
    int nlev;                               // nlev is the number of levels
    double epsilon = 128 * FLT_EPSILON;     // Epsilon used for comparing flaots
    double *y;                              // Solution vector used for spmv
    double *tmp;                            // temp vector for storing old x value, used for sanity checking

    std::vector<int> non_zero;
    util::create_csc(argv[1], argv[2], &Lp, &Li, &Lx, n, &x, non_zero);
    util::create_level_set(n, &Lp, &Li, &jlev, &ilev, nlev, non_zero);

    //Copying right handside vector for sanity check
    tmp = new double[n];              // vector used for sanity check
    for(int i = 0; i < n; i++){
       tmp[i] = x[i]; 
    }

    //TODO: Does a triangular solve
    //sptrsv_csc(n, Lp, Li, Lx, x);
    parallel_sptrsv_csc(n, Lp, Li, Lx, x, nlev, ilev, jlev);
    

    //TODO: sanity check using spmv
    y = new double[n];
    std::fill(y, y+n, 0.);
    //spmv_csc(n, Lp, Li, Lx, x, y);
    parallel_spmv_csc(n, Lp, Li, Lx, x, y);

    
    for(int i = 0; i < n; i++){
        if (!util::AreSame(tmp[i], y[i], epsilon)){
            printf("failed for unknown %d\n", i+1);
            printf("x: %f, y: %f\n", tmp[i], y[i]);
            return 1;
        }
    }
    printf("Passed!\n");

    delete[] Lx;
    delete[] Lp;
    delete[] Li;
    delete[] x;
    delete[] tmp;


    return 0;
}
