#include "matrix_solve.h"

/*      
 * Solves a sparse lower triangular solve (Lx=b) which is stored in CSC format.
 * The sparse matrix is stored in {n, Lp, Li, Lx}
 * x is in/out. It is initially b and at the end it has the solution.
 */

void sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x) {
    int i, j;

    double start_time = omp_get_wtime();
    for (i = 0; i < n; i++) {
        x[i] /= Lx[Lp[i]];
        for (j = Lp[i] + 1; j < Lp[i + 1]; j++) {
            x[Li[j]] -= Lx[j] * x[i];
        }
    }
    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;
    printf("total time: %f\n", elapsed_time);

}

//TODO1: parallel sptrsv_csc
void parallel_sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x, int nlev, int *ilev, int *jlev) {
    int m, j, k, i;

    double start_time = omp_get_wtime();
    for(m = 0; m < nlev; m++){
        #pragma omp parallel for schedule(static)
        for(k = ilev[m]; k < ilev[m+1]; k++){
            i = jlev[k];
            x[i] /= Lx[Lp[i]]; 
            for(j = Lp[i] + 1; j < Lp[i+1]; j++){
                #pragma omp atomic
                x[Li[j]] -= Lx[j] * x[i];
            }
        }
    }
    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;
    printf("total time: %f\n", elapsed_time);

    // Free jlev and ilev since they are no longer needed
    delete[] jlev;
    delete[] ilev;

}

/*
 * Computes y = A*x where A is a sparse matrix {n, Ap, Ai, Ax} and 
 * x and y are vectors
 */
void spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax,
        const double *x, double *y) {

    int i, j;
    for (i = 0; i < n; i++) {
        for (j = Ap[i]; j < Ap[i + 1]; j++) {
            y[Ai[j]] += Ax[j] * x[i];
        }
    }  
}

//TODO2: parallel spmv_csc
void parallel_spmv_csc(int n, const int *Ap, const int *Ai, const double *Ax,
        const double *x, double *y) {

    #pragma omp parallel num_threads(1)
    {
        //double *priv_y = new double[n];
        //std::fill(priv_y, priv_y + n, 0.);

        // Use atomic instead of reduction, reduction seems to pop the 
        // stack sometimes
        #pragma omp parallel for 
        for (int i = 0; i < n; i++) {
            for (int j = Ap[i]; j < Ap[i + 1]; j++) {
                //priv_y[Ai[j]] += Ax[j] * x[i];
                #pragma omp atomic
                y[Ai[j]] += Ax[j] * x[i];
            }
        }
        
        //#pragma omp critical
        //{
        //    for(int i=0;  i< n; ++i) {
        //        y[i] += priv_y[i];
        //    }
        //    delete[] priv_y;
        //}
    }  
}


//TODO3: Serial sptrsv_csr

//TODO4: Serial spmv_csr


