CUDAFLAGS=-O3 -arch=sm_30 -g -ccbin /usr/bin/g++-5 
	
all: CFLAGS= -O3
all: e1

debug: CFLAGS=-g -O3
debug: e1
	gdb e1

e1: matriz.o e1.o complex.o cuda.o
	nvcc $(CUDAFLAGS) -o e1 complex.o e1.o matriz.o cuda.o

e1.o:e1.cpp
	nvcc $(CUDAFLAGS) -dc -c e1.cpp

matriz.o:matriz.cpp matriz.h
	nvcc $(CUDAFLAGS) -dc -x cu -c matriz.cpp

complex.o:complex.cpp complex.h
	nvcc $(CUDAFLAGS) -dc -x cu -c complex.cpp

cuda.o:cuda.cu cuda.cuh
	nvcc $(CUDAFLAGS) -dc -c cuda.cu

clean:
	rm *.o e1
