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
void parallel_sptrsv_csc(int n, int *Lp, int *Li, double *Lx, double *x, std::map<int, std::vector<int>> &lvl_set) {
    int i, j;
    for (i = 0; i < n; i++) {
        x[i] /= Lx[Lp[i]];
        #pragma omp parallel
        {
            #pragma omp parallel for schedule(static) 
            for (j = Lp[i] + 1; j < Lp[i + 1]; j++) {
                    x[Li[j]] -= Lx[j] * x[i];
            }
        }
    }
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

void create_csc(char *matrix, char *b, int *Lp, int *Li, double *Lx, int &n, double *x, std::vector<int> &non_zero, std::map<int, std::vector<int>> &adj_list){
    int R, C, r, c;
    double data;

    // Read answer vector 'b' to find non-zeros 
    std::ifstream bin(b);
    while (bin.peek() == '%') bin.ignore(2048, '\n');  // Ignore comments

    bin >> R >> C >> n;
    x = new double[R]; 
    std::fill(x, x + R, 0.);

    for(int i = 0; i < n; i++){
        bin >> r >> c >> data;
        x[r-1] = data;
        non_zero.push_back(r-1);
    }

    // Read lower triangular matrix
    std::ifstream fin(matrix);
    while (fin.peek() == '%') fin.ignore(2048, '\n');  // Ignore comments

    // Initialize Lx and fill with zero
    fin >> R >> C >> n;   
    Lp = new int[R];                    // Array for storing the diagonal entries' index in Lx
    Li = new int[n];                    // Array for storing the row of entry in Lx
    Lx = new double[n];	                // Array for storing all the non zeros in the matrix

    for(int i = 0; i < n; i++){
        fin >> r >> c >> data;
        Lx[i] = data;
        Li[i] = r-1;
        if (r == c){
          Lp[r-1] = i;
        }
        adj_list[c-1].push_back(r-1);
    }   

    //TODO create working set by iterating over nonzero and Lx and add to a set then iterate over the set instead of R
    /*std::set<int> it_space;
    std::set<int>::iterator it;
    std::pair<std::set<int>::iterator,bool> ret;

    for(int i = 0; i < non_zero.size(); i++){
        for(int j = Lp[non_zero[i]]; j < Lp[non_zero[i]+1]; j++){
            if(it_space.empty()){
                ret = it_space.insert(Li[j]);
            }else{
                ret = it_space.insert(it, Li[j]);
                if (ret.second==false) it=ret.first;
            }   
        }
    }

    int t, l, nlev = 0;
    int *levels = new int[R];
    std::fill(levels, levels+R, 0);

    for(auto i : it_space){
       l = levels[i] + 1; 
       for(int j = Lp[i]; j < Lp[i + 1]; j++){
           t = max(l, levels[Li[j]]);
           nlev = max(nlev, t);
           levels[Li[j]] = t;
       }
    }*/
    
    int t, l, nlev = 0;
    int levels[R] = {0};

    for(int i = 0; i < R; i ++){
       l = levels[i] + 1;     
       for(int j = Lp[i]; j < Lp[i + 1]; j++){
           t = std::max(l, levels[Li[j]]);
           nlev = std::max(nlev, t);
           levels[Li[j]] = t;
       }
    }

    levels[R-1] += 1;
    


    printf("Lx: [");
    for(int i = 0; i < n; i++){
        std::cout << Lx[i] << ", ";
    }
    printf("]\n");

    printf("Lp: [");
    for(int i = 0; i < R; i++){
        std::cout << Lp[i] << ", ";
    }
    printf("]\n");

    printf("Li: [");
    for(int i = 0; i < n; i++){
        std::cout << Li[i] << ", ";
    }
    printf("]\n");

    printf("adj_list: [\n");
    for(int i = 0; i < adj_list.size(); i++){
        std::cout << i+1 << ": ";
        for(int j = 0; j < adj_list[i].size(); j++){
            std::cout << adj_list[i][j] << ", ";
        }
        std::cout << std::endl;
    }
    printf("]\n");

    printf("non_zero [");
    for(int i = 0; i < non_zero.size(); i++){
        std::cout << non_zero[i] << ", ";
    }
    printf("]\n");

    printf("levels: [");
    for(int i = 0; i < R; i++){
        if(levels[i] > 0){
            std::cout << "i[" << i << "]: " << levels[i] << ", ";
        }
    }
    printf("]\n");
}

void create_lvl_set(std:: vector<int> non_zero, std::map<int, std::vector<int>> &adj_list, std::map<int, std::vector<int>> &lvl_set){
    int level[adj_list.size()] = {0};

    #pragma omp parallel
    {
        int priv_lvl[adj_list.size()] = {0};
        std::stack<int> visited;

        // Run dfs on all non zero entry in solution vector 'b'
        #pragma omp for
        for(int i = 0; i < non_zero.size(); i++){
           visited.push(non_zero[i]);

           while(!visited.empty()){
               int curr = visited.top();
               visited.pop();
               priv_lvl[curr]++;

               auto al = adj_list[curr];
               for(auto it : al){
                   if(it != curr){
                       visited.push(it);
                   }
               }
           }
        }

        printf("@@@@@@ OUT AFTER @@@@@@@\n");

        #pragma omp critical
        for (int n = 0; n < adj_list.size(); n++){
            level[n] += priv_lvl[n];
        }
    }

    printf(";alskjfk;asflajfl;s");

    // Scan level array to create lvl set
    for(int i = 0; i < adj_list.size(); i++){
        if(level[i] != 0){
            lvl_set[level[i]-1].push_back(i);
        }
    }

    printf("Level set:\n");
    for(int i = 0; i < lvl_set.size(); i++){
        std::cout << "lvl " << i << ": ";
        for(int j = 0; j < lvl_set[i].size(); j++){
            std::cout << lvl_set[i][j] << " ";
        }
        std::cout << std::endl;
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
    int *Lp, *Li, n;
    std::vector<int> non_zero;
    std::map<int, std::vector<int>> adj_list;
    std::map<int, std::vector<int>> lvl_set;
    create_csc(argv[1], argv[2], Lp, Li, Lx, n, x, non_zero, adj_list);
    printf("@@@@@@@@ done scanning @@@@@@@@\n");
    printf("non_zero sz: %d\n", (int)non_zero.size());
    
    /* 
    printf("non_zero [");
    for(int i = 0; i < non_zero.size(); i++){
        std::cout << non_zero[i] << ", ";
    }
    printf("]\n");
    */
   
   // create_lvl_set(non_zero, adj_list, lvl_set);

 	//TODO: Does a triangular solve
  


 	//TODO: sanity check using spmv

    
 	return 1;
 }
