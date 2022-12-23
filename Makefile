all: openMPI MPI CUDA

openMPI:
	mpic++ -fopenmp -O5 parallelOMPI.cpp -o parallelO.o -lm
MPI:
	mpicc -O5 parallelMPI.cpp -o parallel.o -lm
CUDA:
	bash dvm c prog.c

clean:
	rm parallelO.o parallel.o imperative.o