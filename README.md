This is a working implementation of the lower triangular solve with level sets, but not
yet h-sets. It implements the optimizations described in Parsy paper of reading the RHS
vector to prune the iteration space of unknowns whose values are zero. This is done by
reading the RHS vector and then adding non-zero unknowns to a working set when scanning
the input matrix.


The make file, Makefile, targets most Ubuntu environments. Makefile2 is specifically 
for MacOS since Apple's default clang compiler have trouble compiling openmp. To work
around this I downloaded llvm and have been using the llvm clang on my Mac. 

Instructions:

To run the code do make in the command line. The make file should generate an executable
called 'solve'. Solve takes in two arguements the first one being the name of the matrix
in lower triangular form and second being the RHS vector.

e.g.

$make
$./solve

tests/test.mtx tests/test_b.mtx
