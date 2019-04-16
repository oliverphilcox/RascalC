## MAKEFILE FOR RascalC/grid_power.cpp. This compiles the grid_power.cpp file into the ./power exececutable.

CC = gcc
CFLAGS = -O3 -Wall
CXXFLAGS = -DPOWER -Wall -O3 -DOPENMP
# disable OPENMP to run single threaded
#-DPERIODIC # use this to enable periodic behavior

CXX = g++ -fopenmp -lgomp -std=c++0x -ffast-math

AUNTIE	= power
AOBJS	= grid_power.o 

LD	= g++
LFLAGS	= -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lgomp


main: $(AUNTIE)

$(AUNTIE):	$(AOBJS)
	$(LD) $(AOBJS) $(LFLAGS) -o $(AUNTIE)

clean:
	rm power 
