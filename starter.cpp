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


void print(int* array, int sz){
    printf(" [");
    for(int i = 0; i < sz; i++){
        printf("%d, ", array[i]);
    }
    printf("]\n");
}

/*      
* Solves a sparse lower triangular solve (Lx=b) which is stored in CSC format.
* The sparse matrix is stored in {n, Lp, Li, Lx}
* x is in/out. It is initially b and at the end it has the solution.
*/

void sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x) {
    int i, j;

    //double start_time= omp_get_wtime();
    clock_t begin = clock();
    for (i = 0; i < n; i++) {
        x[i] /= Lx[Lp[i]];
        for (j = Lp[i] + 1; j < Lp[i + 1]; j++) {
            x[Li[j]] -= Lx[j] * x[i];
        }
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("total time: %f\n", time_spent);

}

//TODO1: parallel sptrsv_csc
void parallel_sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double **x, int nlev, int *ilev, int *jlev) {
    int m, j, k, i;

    double start_time= omp_get_wtime();

    for(m = 0; m < nlev; m++){
        #pragma omp parallel for schedule(static)
        for(k = ilev[m]; k < ilev[m+1]; k++){
            i = jlev[k];
            (*x)[i] /= Lx[Lp[i]]; 
            for(j = Lp[i] + 1; j < Lp[i+1]; j++){
                //#pragma omp atomic
                (*x)[Li[j]] -= Lx[j] * (*x)[i];
            }
        }
    }
    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;
    printf("total time: %f\n", elapsed_time);

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


void create_csc(char *matrix, char *b, int **Lp, int **Li, double **Lx, int &n, double **x, std::vector<int> &non_zero){
    int R, C, r, c, num_entry;
    double data;

    // Read answer vector 'b' to find non-zeros 
    std::ifstream bin(b);
    while (bin.peek() == '%') bin.ignore(2048, '\n');  // Ignore comments

    bin >> R >> C >> num_entry;
    *x = new double[R]; 
    int visited[R] = {0};
    std::fill(*x, *x + R, 0.);

    for(int i = 0; i < num_entry; i++){
        bin >> r >> c >> data;
        (*x)[r-1] = data;
        non_zero.push_back(r-1);
        visited[r-1] = 1;
    }

    // Read lower triangular matrix
    std::ifstream fin(matrix);
    while (fin.peek() == '%') fin.ignore(2048, '\n');  // Ignore comments

    // Initialize Lx and fill with zero
    fin >> R >> C >> num_entry;   
    n = R;
    *Lp = new int[n+1];                    // Array for storing the diagonal entries' index in Lx
    *Li = new int[num_entry];                    // Array for storing the row of entry in Lx
    *Lx = new double[num_entry];	                // Array for storing all the non zeros in the matrix
    (*Lp)[n] = num_entry; 

    for(int i = 0; i < num_entry; i++){
        fin >> r >> c >> data;
        (*Lx)[i] = data;
        (*Li)[i] = r-1;
        if (r == c){
          (*Lp)[r-1] = i;
        }

        // Essentially runs dfs adds all affected unknowns to non_zero vector
        if(visited[c-1] && !visited[r-1]){
            visited[r-1] = 1;
            non_zero.push_back(r-1);
        }
    }   

    std::sort(non_zero.begin(), non_zero.end());

}

void create_level_set(int n, int **Lp, int **Li, int **levels, int **jlev, int **ilev, int &nlev, std::vector<int> non_zero){
    int i, t, cnt = 0, l = 0;
    nlev = 0;     

    *levels = new int[n];
    (*jlev) = new int[n];                   // array that denotes the unknowns in ascending order of their levels
    std::fill(*levels, *levels + n, 0);
    std::fill(*jlev, *jlev + n, 0);
    std::vector<int>::iterator it = non_zero.begin();
    for(; it != non_zero.end(); it++){
       i = *it;
       l = (*levels)[i] + 1;     
       for(int j = (*Lp)[i]; j < (*Lp)[i + 1]; j++){
           t = std::max(l, (*levels)[(*Li)[j]]);
           nlev = std::max(nlev, t);
           (*levels)[(*Li)[j]] = t;
       }
    }

    (*levels)[n-1] += 1;
    nlev = std::max(nlev, (*levels)[n-1]);

    (*ilev) = new int[nlev+1];              // array that contains the pointers to the start of the ith level
    std::fill(*ilev, *ilev + nlev + 1, 0);
    (*ilev)[nlev] = n;

    for(int i = 1; i <= nlev; i++){
       (*ilev)[i-1] = cnt; 
       for(int j = 0; j < n; j++){
           if((*levels)[j] == i){
               (*jlev)[cnt] = j;
               cnt++;
           }
       }
    }
}

int main(int argc, char *argv[]){
    
    printf("argc: %d\n", argc);
    if(argc != 3){
        printf("wrong # of args\n");
        return 0;
    }

    //TODO: read a matrix from Suite sparse matrix repo
    double *Lx, *x;
    int *Lp, *Li, *levels, *jlev, *ilev, n, nlev;       // nlev is the number of levels
    std::vector<int> non_zero;
    std::map<int, std::vector<int>> lvl_set;
    create_csc(argv[1], argv[2], &Lp, &Li, &Lx, n, &x, non_zero);
    create_level_set(n, &Lp, &Li, &levels, &jlev, &ilev, nlev, non_zero);

    double y[n] = {0};
    for(int i = 0; i < n; i++){
       y[i] = x[i]; 
    }

    //TODO: Does a triangular solve
    //sptrsv_csc(n, Lp, Li, Lx, x);
    parallel_sptrsv_csc(n, Lp, Li, Lx, &x, nlev, ilev, jlev);

    //TODO: sanity check using spmv
    double tmp[n] = {0};
    spmv_csc(n, Lp, Li, Lx, x, tmp);
    
    printf("\n");
    for(int i = 0; i < n; i++){
        if (tmp[i] != y[i]){
            printf("\n i: %d fuck\n", i);
            return 0;
        }
    }

    return 1;
}
