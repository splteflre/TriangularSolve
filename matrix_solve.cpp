#include "starter.h"

/*      
 * Solves a sparse lower triangular solve (Lx=b) which is stored in CSC format.
 * The sparse matrix is stored in {n, Lp, Li, Lx}
 * x is in/out. It is initially b and at the end it has the solution.
 */

void sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double **x) {
    int i, j;

    double start_time= omp_get_wtime();
    for (i = 0; i < n; i++) {
        (*x)[i] /= Lx[Lp[i]];
        for (j = Lp[i] + 1; j < Lp[i + 1]; j++) {
            (*x)[Li[j]] -= Lx[j] * (*x)[i];
        }
    }
    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;
    printf("total time: %f\n", elapsed_time);

}

//TODO1: parallel sptrsv_csc
void parallel_sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x, int nlev, int *ilev, int *jlev) {
    int m, j, k, i;

    double start_time= omp_get_wtime();
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


//TODO3: Serial sptrsv_csr

//TODO4: Serial spmv_csr

/*
 * Parse and read the matrix files described by the variable 'matrix' and 'b'
 */
void create_csc(char *matrix, char *b, int **Lp, int **Li, double **Lx, int &n, double **x, std::vector<int> &non_zero){
    int R, C, r, c, num_entry;
    double data;
    int *visited;

    // Read answer vector 'b' to find non-zeros 
    std::ifstream bin(b);
    while (bin.peek() == '%') bin.ignore(2048, '\n');  // Ignore comments

    bin >> R >> C >> num_entry;
    *x = new double[R]; 
    visited = new int[R];
    std::fill(*x, *x + R, 0.);
    std::fill(visited, visited + R, 0);

    for(int i = 0; i < num_entry; i++){
        bin >> r >> c >> data;
        (*x)[r-1] = data;
        non_zero.push_back(r-1);
        visited[r-1] = 1;
    }
    bin.close();

    // Read lower triangular matrix
    std::ifstream fin(matrix);
    while (fin.peek() == '%') fin.ignore(2048, '\n');  // Ignore comments

    // Initialize Lx and fill with zero
    fin >> R >> C >> num_entry;   
    n = R;
    *Lp = new int[n+1];                    
    *Li = new int[num_entry];                
    *Lx = new double[num_entry];	        
    (*Lp)[n] = num_entry; 

    for(int i = 0; i < num_entry; i++){
        fin >> r >> c >> data;
        (*Lx)[i] = data;
        (*Li)[i] = r-1;
        if (r == c){
            (*Lp)[r-1] = i;
        }

        // Filter out all unknowns who's values are zero, faster to do this than run dfs
        if(visited[c-1] && !visited[r-1]){
            visited[r-1] = 1;
            non_zero.push_back(r-1);
        }
    }   
    fin.close();

    std::sort(non_zero.begin(), non_zero.end());

}

/*
 * Create level sets as described in the sympiler paper
 *
 */
void create_level_set(int n, int **Lp, int **Li, int **jlev, int **ilev, int &nlev, std::vector<int> non_zero){
    int i, j, t, cnt = 0, l = 0;
    int *levels = new int[n];               // Array that maps the known at index i to its level 
    nlev = 0;     

    (*jlev) = new int[n];
    std::fill(levels, levels + n, 0);
    std::fill(*jlev, *jlev + n, 0);

    // Create the level set for all unknowns whose position in the right handside vector is either
    // a non zero or is affected by a non zero unknown. tldr: create level set for all non zero unknown
    for(auto it = non_zero.begin(); it != non_zero.end(); it++){
        i = *it;
        l = levels[i] + 1;     
        for(int j = (*Lp)[i]; j < (*Lp)[i + 1]; j++){
            t = std::max(l, levels[(*Li)[j]]);
            nlev = std::max(nlev, t);
            levels[(*Li)[j]] = t;
        }
    }
    levels[n-1] += 1;
    nlev = std::max(nlev, levels[n-1]);

    (*ilev) = new int[nlev+1];
    std::fill(*ilev, *ilev + nlev + 1, 0);
    (*ilev)[nlev] = n;

    // Create the level set pointers and array of unknowns in ascending order of their levels
    for(i = 1; i <= nlev; i++){
        (*ilev)[i-1] = cnt; 
        for(j = 0; j < n; j++){
            if(levels[j] == i){
                (*jlev)[cnt] = j;
                cnt++;
            }
        }
    }
    delete[] levels;
}

