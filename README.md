# CSC418-bonus


Overview:

The two papers describe a specific technique to improve sparse matrix computations, specifically, solving lower triangular matrices with forward substitution in parallel. One paper goes in depth about heuristics to improve parallelism for when matrices are large, while the other paper talks about code generation for specific matrices. This is actually really cool! Although matrices may be different numerically, the sparsity of the matrices tend to stay the same. Basically, the numbers themselves may be different but where the numbers are located in the actual matrix tends to stay the same. This allows for code targeted specific matrices to be used repeatedly even if the numbers themselves are different.


The technique described in the paper first creates a DAG(directed acyclic graph) where each vertice is an unknown and each edge describes if the unknown at the head of the edge requires the unknown at the tail of the edge to solve it. For example, if X1 -> X2 this means that we need to solve X1 before solving X2.  We then analyze the right hand side vector, aka the ‘b’ in Lx = b, to see where the non zeros are located. We use those non zeros to filter out unknowns that aren't needed and to create “level sets” (this is for running it in parallel). Level sets tell us which columns we can solve in parallel.


I implemented both the serial and parallel version of solving the lower triangular matrices for both CSC and CSR storage formats. I also implemented a way to check if it worked in serial and parallel for CSC but only serial for CSR. I did both CSC and CSR formats since I implemented CSC first but it turned out to be super slow and I couldn’t get it to go any faster so I did the CSR to see if it would be better.


Pitfalls:

The lock used in the CSC implementation KILLS parallelism. To put it in perspective for a 503625 by 503625 matrix with 9027150 non zero entries my parallel ran for around ~20 seconds but the serial one took about 1 second. My CSR did improve however I think there’s still a lot of false sharing which kills some of the parallelism.


