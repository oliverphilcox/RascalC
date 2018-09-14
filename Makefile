#CXXFLAGS = -O2 -Wall 
#CXXFLAGS = -g -Wal

#CC = icc -ipo
#CXX = icc -liomp5 -openmp

CC = gcc
CFLAGS = -O3 -Wall
CXXFLAGS = -O3 -Wall -DOPENMP -DPERIODIC
#CXXFLAGS = -Wall -DOPENMP -O3 -unroll-aggressive -opt-prefetch -ipo

CXX = g++ -fopenmp -lgomp -std=c++0x -ffast-math
#CXXFLAGS = -O2 -DOPENMP #-openmp

AUNTIE	= cov
AOBJS	= grid_covariance.o ./cubature/hcubature.o ./ransampl/ransampl.o

LD	= g++
LFLAGS	= -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lgomp


main: $(AUNTIE)

$(AUNTIE):	$(AOBJS)
	$(LD) $(AOBJS) $(LFLAGS) -o $(AUNTIE)

#grid_covariance.o: ../threept/interp2d/interp2d_spline.h

#default: grid_covariance

clean:
	rm grid_covariance
