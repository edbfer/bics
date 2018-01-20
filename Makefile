CUDAFLAGS=-O3 -arch=sm_30 -g -ccbin /usr/bin/g++-5 
CFLAGS=-g -O3
all: CFLAGS= -O3
all: e1

e4: e4

debug: CFLAGS=-g -O3
debug: e1
	gdb e1

e4: e4.o
	g++ $(CFLAGS) -o e4 e4.o

e4.o: g6s6c1a.cpp
	g++ $(CFLAGS) -o e4.o -c g6s6c1a.cpp

e1: e1.o
	g++ $(CUDAFLAGS) -o e1 e1.o

e1.o:e1.cpp
	g++ $(CFLAGS) -c e1.cpp


clean:
	rm *.o e1
