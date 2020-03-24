This is a working implementation of the lower triangular solve with level sets, but not
yet h-sets. It implements the optimizations described in Parsy paper of reading the RHS
vector to prune the iteration space of unknowns whose values are zero. This is done by
reading the RHS vector and then adding non-zero unknowns to a working set when scanning
the input matrix.


The make file, Makefile, targets most Ubuntu environments. Makefile2 is specifically 
for MacOS since Apple's default clang compiler have trouble compiling openmp. To work
around this I downloaded llvm and have been using the llvm clang on my Mac. 

