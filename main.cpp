#include "util/util.h"
#include "matrix_solve/matrix_solve.h"


int main(int argc, char *argv[]){
    
    printf("argc: %d\n", argc);
    if(argc != 4){
        printf("wrong # of args\n");
        printf("Usage: solve <Lower Triangular Matrix> <RHS vector> <Matrix format CSC/CSR>\n");
        return 0;
    }

    //TODO: read a matrix from Suite sparse matrix repo
    double *Lx; 
    int *Lp, *Li;                           // Lx, Lp, and Li arrays as described in the handout
    double *x;                              // Right hand side vector, b
    double *d;                              // Array for diagonal entries in L
    int n;                                  // Dimension of the matrix read from input file
    int *jlev;                              // Array that lists the unknowns in ascending order of their levels
    int *ilev;                              // Array that contains the pointers to the levels in jlev
    int nlev;                               // nlev is the number of levels
    double epsilon = 128 * FLT_EPSILON;     // Epsilon used for comparing flaots
    double *tmp_x;                           // Temp vector for storing old x value, used for sanity checking
    double *y;                               // Solution vector for spmv
    std::vector<int> non_zero;              // Vector for non zero entries in b answer vector
    char* format = argv[3];                 // Storage format, either CSC or CSR

    util::init_arrays(argv[1], argv[2], &Lp, &Li, &Lx, n, &x, non_zero, format, &d);
    util::create_level_set(n, &Lp, &Li, &jlev, &ilev, nlev, non_zero, format);

    //Copying right handside vector for sanity check
    tmp_x = new double[n];                    // vector used for sanity check
    for(int i = 0; i < n; i++){
       tmp_x[i] = x[i]; 
    }

    y = new double[n];
    std::fill(y, y+n, 0.);

    if(strcmp(format, "CSR") == 0){
        // Triangle solve with CSR
        //sptrsv_csr(n, Lp, Li, Lx, x, d);
        parallel_sptrsv_csr(n, Lp, Li, Lx, tmp_x, nlev, ilev, jlev, d);

        // Sanity check using spmv
        spmv_csr(n, Lp, Li, Lx, tmp_x, y);
    }else{
        // Triangle solve with CSC
        //sptrsv_csc(n, Lp, Li, Lx, x);
        parallel_sptrsv_csc(n, Lp, Li, Lx, tmp_x, nlev, ilev, jlev);

        // Sanity check using spmv
        spmv_csc(n, Lp, Li, Lx, tmp_x, y);
        //parallel_spmv_csc(n, Lp, Li, Lx, tmp_x, y);
    }
    
    
    for(int i = 0; i < n; i++){
        if (!util::AreSame(x[i], y[i], epsilon)){
            printf("failed for unknown %d\n", i+1);
            printf("x: %f, y: %f\n", x[i], y[i]);
            return 1;
        }
    }
    printf("Passed!\n");
    
    delete[] Lx;
    delete[] Lp;
    delete[] Li;
    delete[] tmp_x;
    delete[] x;
    delete[] y;
    if(strcmp(format, "CSR") == 0)
        delete[] d;
    return 0;
}
