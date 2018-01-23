CUDAFLAGS=-O3 -arch=sm_30 -g -ccbin /usr/bin/g++-5 
CFLAGS=-g

all: e5

e4: e4

debug: CFLAGS=-g -O3
debug: e1
	gdb e1

helper: helper.c
	g++ $(CFLAGS) -o helper helper.c

e4: e4.o
	g++ $(CFLAGS) -o e4 e4.o

e4.o: g6s6c1a.cpp
	g++ $(CFLAGS) -o e4.o -c g6s6c1a.cpp

e5: e5.o
	g++ $(CFLAGS) -o e5 e5.o

e5.o: e5.cpp
	g++ $(CFLAGS) -o e5.o -c e5.cpp

e1: e1.o
	g++ $(CFLAGS) -o e1 e1.o

e1.o:e1.cpp
	g++ $(CFLAGS) -c e1.cpp


clean:
	rm *.o e1 e5