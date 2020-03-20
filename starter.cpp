#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm> 
#include <vector>
#include <stack>
#include <set>
#include <map>
#include <omp.h>


void print(double* array, int sz){
    printf(" [");
    for(int i = 0; i < sz; i++){
        printf("%f, ", array[i]);
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
  for (i = 0; i < n; i++) {
   x[i] /= Lx[Lp[i]];
   for (j = Lp[i] + 1; j < Lp[i + 1]; j++) {
    x[Li[j]] -= Lx[j] * x[i];
   }
  }
 }

//TODO1: parallel sptrsv_csc
void parallel_sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x, int nlev, int *ilev, int *jlev) {
    int m, j, k, i;

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
    printf("x: ");
    print(x, n);
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


void create_csc(char *matrix, char *b, int **Lp, int **Li, double **Lx, int &n, double **x, int **levels, int **jlev, int **ilev, int &nlev){
    int R, C, r, c, num_entry;
    double data;

    // Read answer vector 'b' to find non-zeros 
    std::ifstream bin(b);
    while (bin.peek() == '%') bin.ignore(2048, '\n');  // Ignore comments

    bin >> R >> C >> num_entry;
    *x = new double[R]; 
    std::fill(*x, *x + R, 0.);

    for(int i = 0; i < num_entry; i++){
        bin >> r >> c >> data;
        (*x)[r-1] = data;
        //non_zero.push_back(r-1);
    }

    // Read lower triangular matrix
    std::ifstream fin(matrix);
    while (fin.peek() == '%') fin.ignore(2048, '\n');  // Ignore comments

    // Initialize Lx and fill with zero
    fin >> R >> C >> num_entry;   
    n = R;
    *Lp = new int[n];                    // Array for storing the diagonal entries' index in Lx
    *Li = new int[num_entry];                    // Array for storing the row of entry in Lx
    *Lx = new double[num_entry];	                // Array for storing all the non zeros in the matrix

    for(int i = 0; i < num_entry; i++){
        fin >> r >> c >> data;
        (*Lx)[i] = data;
        (*Li)[i] = r-1;
        if (r == c){
          (*Lp)[r-1] = i;
        }
    }   

    *levels = new int[n];
    (*jlev) = new int[n];       // array that denotes the unknowns in ascending order of their levels
    (*ilev) = new int[nlev+1];  // array that contains the pointers to the start of the ith level
    std::fill(*jlev, *jlev+n, 0);
    std::fill(*ilev, *ilev+nlev+1, 0);
}

void create_level_set(int n, int **Lp, int **Li, int **levels, int **jlev, int **ilev, int &nlev){
    int t, cnt = 0, l = 0;
    nlev = 0;     

    for(int i = 0; i < n; i++){
       l = (*levels)[i] + 1;     
       for(int j = (*Lp)[i]; j < (*Lp)[i + 1]; j++){
           t = std::max(l, (*levels)[(*Li)[j]]);
           nlev = std::max(nlev, t);
           (*levels)[(*Li)[j]] = t;
       }
    }
    (*levels)[n-1] += 1;
    nlev = std::max(nlev, (*levels)[n-1]);
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
    create_csc(argv[1], argv[2], &Lp, &Li, &Lx, n, &x, &levels, &jlev, &ilev, nlev);
    create_level_set(n, &Lp, &Li, &levels, &jlev, &ilev, nlev);
    parallel_sptrsv_csc(n, Lp, Li, Lx, x, nlev, ilev, jlev);
    
    //parallel_sptrsv_csc(n, Lp, Li, Lx, x, lvl_set);


    //TODO: Does a triangular solve
    //parallel_sptrsv_csc(n, Lp, Li, Lx, x, lvl_set); 



    //TODO: sanity check using spmv


    return 1;
}
