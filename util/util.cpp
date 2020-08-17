#include "util.h"

namespace util
{

    /*
     * Default double comparison is not accurate enough using the code comparison found here instead:
     * https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
     */
    bool AreSame(double a, double b, double epsilon)
    {
        // Testing different floating point comparison methods
        //return a == b;
        //double relth = FLT_MIN;
        //double diff = std::abs(a-b);
        //double norm = std::min((abs(a) + abs(b)), std::numeric_limits<double>::max());
        //return diff < std::max(relth, epsilon * norm);

        return fabs(a - b) < epsilon;
    }

    /*
     * Parse and read the matrix files described by the variable 'matrix' and 'b'.
     * Depending on 'format' either use the row or column for Li
     */
    void init_arrays(char *matrix, char *b, int **Lp, int **Li, double **Lx, int &n, double **x, std::vector<int> &non_zero, char *format, double **d)
    {
        int R, C, r, c, num_entry;
        double data;
        int *visited;

        // Read answer vector 'b' to find non-zeros 
        std::ifstream bin(b);
        while (bin.peek() == '%') bin.ignore(2048, '\n');  // Ignore comments

        // Read 'b' vector and initialize arrays
        bin >> R >> C >> num_entry;
        *x = new double[R]; 
        visited = new int[R];
        std::fill(*x, *x + R, 0.);
        std::fill(visited, visited + R, 0);

        for(int i = 0; i < num_entry; i++){
            bin >> r >> c >> data;
            (*x)[r-1] = data;
            visited[r-1] = 1;        // Add non-zero RHS entry to visited array, used for prunning iteration space
            non_zero.push_back(r-1);
        }
        bin.close();

        // Read lower triangular matrix
        std::ifstream fin(matrix);
        while (fin.peek() == '%') fin.ignore(2048, '\n');  // Ignore comments

        // Initialize arrays used for csc and fill with zero
        fin >> R >> C >> num_entry;   
        n = R;
        *Lp = new int[n+1];                    
        *Li = new int[num_entry];                
        *Lx = new double[num_entry];	        
        std::fill(*Lp, *Lp + n + 1, 0);
        std::fill(*Li, *Li + num_entry, 0);
        std::fill(*Lx, *Lx + num_entry, 0.0);
        (*Lp)[n] = num_entry; 

        // This might seem dumb but I think its worth since 'n' can get REALLY big
        if(strcmp(format, "CSC") == 0){
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
        }else{
            int entry_idx = 0;
            int diag_idx  = 0;
            bool first_entry = true;
            *d = new double[n];
            std::fill(*d, *d + n, 0.0);


            for(int i = 0; i < num_entry; i++){
                fin >> r >> c >> data;
                (*Lx)[i] = data;
                (*Li)[i] = c-1;
                
                if(first_entry){
                    first_entry = false;
                    (*Lp)[entry_idx] = i;
                    entry_idx += 1;
                }

                if(r == c){
                    first_entry = true;
                    (*d)[diag_idx] = data;
                    diag_idx += 1;
                }

                // Filter out all unknowns who's values are zero, faster to do this than run dfs
                if(visited[c-1] && !visited[r-1]){
                    visited[r-1] = 1;
                    non_zero.push_back(r-1);
                }
            }   
        }
        fin.close();

        std::sort(non_zero.begin(), non_zero.end());

    }

    /*
     * Create level sets as described in the sympiler paper
     *
     */
    void create_level_set(int n, int **Lp, int **Li, int **jlev, int **ilev, int &nlev, std::vector<int> non_zero, char* format)
    {
        int i, j, t, cnt = 0, l = 0;
        int *levels = new int[n];               // Array that maps the unknown at index i to its level 
        nlev = 0;     

        (*jlev) = new int[n];
        std::fill(levels, levels + n, 0);
        std::fill(*jlev, *jlev + n, 0);

        // Create the level set for all unknowns whose position in the right handside vector is either
        // a non zero or is affected by a non zero unknown. tldr: create level set for all non zero unknown
        
        if(strcmp(format, "CSC") == 0){
            std::vector<int>::iterator it = non_zero.begin();
            for(; it != non_zero.end(); it++){
                i = *it;
                l = levels[i] + 1;     
                for(j = (*Lp)[i]; j < (*Lp)[i + 1]; j++){
                    t = std::max(l, levels[(*Li)[j]]);
                    levels[(*Li)[j]] = t;
                    nlev = std::max(nlev, t);
                }
            }
            //levels[n-1] += 1;   
            nlev = std::max(nlev, levels[n-1]);
        }else{
            std::vector<int>::iterator it = non_zero.begin();
            for(; it != non_zero.end(); it++){
                i = *it;
                l = -1;
                for(j = (*Lp)[i]; j < (*Lp)[i + 1]; j++){
                    l = std::max(l, levels[(*Li)[j]]+1);
                }
                levels[i] = l;
                nlev = std::max(nlev, l);
            }
            nlev = std::max(nlev, levels[n-1]);
        }

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

};
