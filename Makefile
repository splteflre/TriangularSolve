# Makefile for Lower Triangular solve

# *****************************************************
# Variables to control Makefile operation

CPP = g++
CPPFLAGS = -Wall -O3 -fopenmp

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
	rm -f main.o matrix_solve.o util.o solve
