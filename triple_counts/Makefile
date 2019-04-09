## MAKEFILE FOR RascalC. This compiles the triple_counts.cpp file into the ./triple exececutable.

CC = gcc
CFLAGS = -Wall -O3 -Wall
CXXFLAGS = -Wall -O3 -DOPENMP
# disable OPENMP to run single threaded
#-DPERIODIC # use this to enable periodic behavior

CXX = g++ -fopenmp -lgomp -std=c++0x -ffast-math

AUNTIE	= triple
AOBJS	= triple_counts.o ../cubature/hcubature.o ../ransampl/ransampl.o

LD	= g++
LFLAGS	= -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lgomp


main: $(AUNTIE)

$(AUNTIE):	$(AOBJS)
	$(LD) $(AOBJS) $(LFLAGS) -o $(AUNTIE)

clean:
	rm triple 
