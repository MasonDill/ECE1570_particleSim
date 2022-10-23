#
CC = g++ -std=gnu++0x
MPCC = mpic++
#Note: this is the flag for gnu compilers. Change this to -openmp for Intel compilers
OPENMP = -fopenmp
MPI = 
CFLAGS = 
LIBS = -lm


TARGETS = serial openmp autograder qtree_omp

all:	$(TARGETS)

serial: serial.o common.o quadtree.o
	$(CC) -o $@ serial.o common.o quadtree.o -lm -g
autograder: autograder.o common.o
	$(CC) -o $@ autograder.o common.o -lm -g
openmp: openmp.o common.o
	$(CC) -o $@ $(OPENMP) openmp.o common.o -lm -g
mpi: mpi.o common.o
	$(MPCC) -o $@ $(MPILIBS) mpi.o common.o -lm -g
qtree_omp: qtree_omp.o common.o quadtree.o
	$(CC) -o $@ $(OPENMP) qtree_omp.o common.o quadtree.o -lm -g

autograder.o: autograder.cpp common.h
	$(CC) -c $(CFLAGS) autograder.cpp -g
openmp.o: openmp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) openmp.cpp -g
serial.o: serial.cpp common.h
	$(CC) -c $(CFLAGS) serial.cpp -g
mpi.o: mpi.cpp common.h
	$(MPCC) -c $(MPI) $(CFLAGS) mpi.cpp -g
common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp -g
quadtree.o: quadtree.cpp quadtree.hpp
	$(CC) -c $(CFLAGS) quadtree.cpp -g
qtree_omp.o: qtree_omp.cpp common.h
	$(CC) -c $(OPENMP) $(CFLAGS) qtree_omp.cpp -g

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt
