# Makefile for Lower Triangular solve

# *****************************************************
# Variables to control Makefile operation
CPP = /usr/local/opt/llvm/bin/clang++
CPPFLAGS = -I/usr/local/opt/llvm/include -fopenmp
LDFLAGS = -L/usr/local/opt/llvm/lib

# ****************************************************
solve: main.o matrix_solve.o util.o
	$(CPP) $(CPPFLAGS) -o solve main.o matrix_solve.o util.o 

main.o: main.cpp util/util.h matrix_solve/matrix_solve.h
	$(CPP) $(CPPFLAGS) -c main.cpp 

matrix_solve.o: matrix_solve/matrix_solve.cpp matrix_solve/matrix_solve.h
	$(CPP) $(CPPFLAGS) -c matrix_solve/matrix_solve.cpp

util.o: util/util.cpp util/util.h
	$(CPP) $(CPPFLAGS) -c util/util.cpp

clean:
	rm main.o matrix_solve.o util.o solve
