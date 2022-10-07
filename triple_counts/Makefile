## MAKEFILE FOR RascalC. This compiles the triple_counts.cpp file into the ./triple exececutable.

CC = gcc
CFLAGS = -Wall -O3 -MMD
CXXFLAGS = -Wall -O3 -MMD -DOPENMP
# DOPENMP # use this to run multi threaded with OpenMP
#-DPERIODIC # use this to enable periodic behavior

CXX = g++ -fopenmp -lgomp -std=c++0x -ffast-math $(shell pkg-config --cflags gsl)

AUNTIE	= triple
AOBJS	= triple_counts.o ../cubature/hcubature.o ../ransampl/ransampl.o
ADEPS   = ${AOBJS:.o=.d}

LD	= g++
LFLAGS	= -L/usr/local/lib -L/usr/lib/x86_64-linux-gnu $(shell pkg-config --libs gsl) -lgomp

.PHONY: main clean

main: $(AUNTIE)

$(AUNTIE):	$(AOBJS) Makefile
	$(LD) $(AOBJS) $(LFLAGS) -o $(AUNTIE)

clean:
	rm -f ${AUNTIE} ${AOBJS} ${ADEPS}

$(AOBJS): Makefile
-include ${ADEPS}
