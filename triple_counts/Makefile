## MAKEFILE FOR RascalC. This compiles the triple_counts.cpp file into the ./triple exececutable.

CC = gcc
CFLAGS = -O3 -Wall -MMD
CXX = g++
CXXFLAGS	= -O3 -Wall -MMD -std=c++11 -ffast-math $(shell pkg-config --cflags gsl)
CXXFLAGS	+= -DOPENMP -DTRIPLE -DPRINTPERCENTS
#-DOPENMP  # use this to run multi-threaded with OPENMP
#-DPERIODIC # use this to enable periodic behavior
#-DTRIPLE # need this to have mumin=-1 and not 0 to be consistent with the 3PCF covariance
#-DPRINTPERCENTS # use this to print percentage of progress in each loop. This can be a lot of output

LFLAGS	= $(shell pkg-config --libs gsl) # common part

# Known OS-specific choices
ifeq ($(shell uname -s),Darwin)
# Here we load the Mac OpenMP. Tested after installation commands:
# brew install libomp
# This may need to be modified with a different installation
ifndef HOMEBREW_PREFIX
HOMEBREW_PREFIX = /usr/local
endif
CXXFLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
LFLAGS	+= -L$(HOMEBREW_PREFIX)/opt/libomp/lib -lomp
else
# default (Linux) case
CXXFLAGS += -fopenmp
LFLAGS	+= -lgomp
endif

LD = $(CXX)

AUNTIE	= triple
AOBJS	= triple_counts.o ../cubature/hcubature.o ../ransampl/ransampl.o
ADEPS   = ${AOBJS:.o=.d}

.PHONY: main clean

main: $(AUNTIE)

$(AUNTIE):	$(AOBJS) Makefile
	$(LD) $(AOBJS) $(LFLAGS) -o $(AUNTIE)

clean:
	rm -f ${AUNTIE} ${AOBJS} ${ADEPS}

$(AOBJS): Makefile
-include ${ADEPS}
